"""Training utilities for PDU autoencoders.

This module provides high-level training functions for autoencoders, handling:
- Data loading and batching (single-process and distributed)
- Optimization and loss computation
- GPU/CPU device management
- Distributed Data Parallel (DDP) training across multiple CPUs/GPUs
- Early stopping and checkpointing
- Training logging and statistics

Example (single-process):
    >>> import numpy as np
    >>> from pdusearch.models import train_autoencoder, TrainingConfig
    >>>
    >>> X = np.random.randn(10000, 900).astype(np.float32)
    >>> config = TrainingConfig(latent_dim=16, epochs=10, batch_size=32)
    >>> model, Z, losses = train_autoencoder(X, config)
    >>> print(f"Final loss: {losses[-1]:.4f}")

Example (DDP on ARCHER2):
    >>> # Run with torchrun:
    >>> # torchrun --nproc_per_node=128 bin/train_ddp.py --aa L --encoding aa_20
"""

from typing import Tuple, Optional, Dict, List
import numpy as np
import logging
import os
import torch
import torch.nn as nn
import torch.distributed as dist
from torch.utils.data import DataLoader, TensorDataset, DistributedSampler

from pdusearch.config import TrainingConfig
from .autoencoder import Autoencoder, VariationalAutoencoder

logger = logging.getLogger(__name__)

# Detect DDP environment
def _is_distributed() -> bool:
    """Check if running under distributed training."""
    return "RANK" in os.environ and "WORLD_SIZE" in os.environ

def _get_rank() -> int:
    """Get current process rank."""
    return int(os.environ.get("RANK", 0))

def _get_world_size() -> int:
    """Get total number of processes."""
    return int(os.environ.get("WORLD_SIZE", 1))


def train_autoencoder(
    X: np.ndarray,
    config: TrainingConfig,
    model_type: str = "standard",
    validate_x: Optional[np.ndarray] = None,
    verbose: bool = True,
) -> Tuple[nn.Module, np.ndarray, List[float]]:
    """Train autoencoder on feature data.

    Comprehensive training function handling device management, optimization,
    loss tracking, and embedding extraction. Returns trained model, learned
    embeddings, and training history.

    Args:
        X: Training feature matrix of shape (n_samples, input_dim) as float32
        config: TrainingConfig instance with hyperparameters
        model_type: Type of autoencoder ('standard' or 'variational')
        validate_x: Optional validation data for early stopping
        verbose: Whether to log training progress

    Returns:
        Tuple of (model, Z, losses) where:
            - model: Trained nn.Module (on CPU)
            - Z: Learned embeddings of shape (n_samples, latent_dim) as float32
            - losses: List of loss values per epoch

    Raises:
        ValueError: If input data invalid or config inconsistent

    Example:
        >>> X = np.random.randn(10000, 900).astype(np.float32)
        >>> config = TrainingConfig(latent_dim=16, epochs=10, batch_size=32)
        >>> model, Z, losses = train_autoencoder(X, config)
        >>> print(Z.shape)  # (10000, 16)
    """
    # Validate input
    if X.ndim != 2:
        raise ValueError(f"X must be 2D, got shape {X.shape}")
    if X.dtype != np.float32:
        X = X.astype(np.float32)
    if np.any(np.isnan(X)) or np.any(np.isinf(X)):
        raise ValueError("X contains NaN or inf values")

    n_samples, input_dim = X.shape

    if verbose:
        logger.info(
            f"Training {model_type} autoencoder: "
            f"{input_dim}D → {config.latent_dim}D "
            f"({n_samples:,} samples, {config.epochs} epochs)"
        )

    # Get device
    device = torch.device(config.device)
    if verbose:
        logger.info(f"  Device: {device}")

    # Create model
    if model_type == "standard":
        model = Autoencoder(
            input_dim=input_dim,
            latent_dim=config.latent_dim,
        )
    elif model_type == "variational":
        model = VariationalAutoencoder(
            input_dim=input_dim,
            latent_dim=config.latent_dim,
            beta=0.1,
        )
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    model = model.to(device)

    if verbose:
        logger.info(f"  Model: {model.count_parameters():,} parameters")

    # Create data loader
    X_tensor = torch.FloatTensor(X).to(device)
    dataset = TensorDataset(X_tensor)
    loader = DataLoader(
        dataset,
        batch_size=config.batch_size,
        shuffle=True,
        pin_memory=(device.type == "cuda"),
    )

    # Optimizer
    if config.optimizer.lower() == "adam":
        optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)
    elif config.optimizer.lower() == "sgd":
        optimizer = torch.optim.SGD(model.parameters(), lr=config.learning_rate)
    elif config.optimizer.lower() == "rmsprop":
        optimizer = torch.optim.RMSprop(model.parameters(), lr=config.learning_rate)
    else:
        raise ValueError(f"Unknown optimizer: {config.optimizer}")

    # Loss function
    if config.loss_function.lower() == "mse":
        loss_fn = nn.MSELoss(reduction="mean")
    elif config.loss_function.lower() == "bce":
        loss_fn = nn.BCELoss(reduction="mean")
    else:
        raise ValueError(f"Unknown loss: {config.loss_function}")

    if verbose:
        logger.info(f"  Optimizer: {config.optimizer}, lr={config.learning_rate}")
        logger.info(f"  Loss: {config.loss_function}")

    # Training loop
    losses = []
    best_val_loss = float("inf")
    patience_counter = 0

    for epoch in range(config.epochs):
        model.train()
        epoch_losses = []

        for (batch_X,) in loader:
            optimizer.zero_grad()

            if model_type == "standard":
                batch_recon = model(batch_X)
                loss = loss_fn(batch_recon, batch_X)
            else:  # variational
                batch_recon, mean, log_var = model(batch_X)
                recon_loss = loss_fn(batch_recon, batch_X)
                # KL divergence: -0.5 * sum(1 + log(σ²) - μ² - σ²)
                kl_loss = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
                kl_loss = kl_loss / batch_X.shape[0]  # Normalize by batch size
                loss = recon_loss + model.beta * kl_loss

            loss.backward()
            optimizer.step()
            epoch_losses.append(loss.item())

        epoch_loss = np.mean(epoch_losses)
        losses.append(epoch_loss)

        # Validation check for early stopping
        if validate_x is not None and config.early_stopping_patience is not None:
            model.eval()
            with torch.no_grad():
                validate_x_tensor = torch.FloatTensor(validate_x).to(device)
                if model_type == "standard":
                    val_recon = model(validate_x_tensor)
                    val_loss = loss_fn(val_recon, validate_x_tensor).item()
                else:
                    val_recon, val_mean, val_log_var = model(validate_x_tensor)
                    val_recon_loss = loss_fn(val_recon, validate_x_tensor)
                    val_kl_loss = -0.5 * torch.sum(1 + val_log_var - val_mean.pow(2) - val_log_var.exp())
                    val_kl_loss = val_kl_loss / validate_x_tensor.shape[0]
                    val_loss = (val_recon_loss + model.beta * val_kl_loss).item()

                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    patience_counter = 0
                else:
                    patience_counter += 1

                if patience_counter >= config.early_stopping_patience:
                    if verbose:
                        logger.info(
                            f"Early stopping at epoch {epoch + 1}/"
                            f"{config.epochs} (val_loss={val_loss:.4f})"
                        )
                    break

        # Periodic logging
        if verbose and (epoch + 1) % max(1, config.epochs // 5) == 0:
            logger.info(f"  Epoch {epoch + 1:3d}/{config.epochs}: loss={epoch_loss:.4f}")

    if verbose:
        logger.info(f"✓ Training complete: final loss={losses[-1]:.4f}")

    # Extract embeddings
    model.eval()
    with torch.no_grad():
        if model_type == "standard":
            Z_tensor = model.encode(X_tensor)
        else:
            mean, _ = model.encode(X_tensor)
            Z_tensor = mean  # Use mean for embeddings (no sampling in eval)

    Z = Z_tensor.cpu().numpy().astype(np.float32)

    # Move model to CPU for return
    model = model.cpu()

    return model, Z, losses


def train_autoencoder_ddp(
    X: np.ndarray,
    config: TrainingConfig,
    model_type: str = "standard",
    validate_x: Optional[np.ndarray] = None,
) -> Tuple[nn.Module, np.ndarray, List[float]]:
    """Train autoencoder using Distributed Data Parallel (DDP).

    This function is designed to run under torchrun with multiple processes.
    Automatically detects DDP environment via RANK/WORLD_SIZE env vars.

    Args:
        X: Training feature matrix of shape (n_samples, input_dim) as float32
        config: TrainingConfig instance with hyperparameters
        model_type: Type of autoencoder ('standard' or 'variational')
        validate_x: Optional validation data for early stopping

    Returns:
        Tuple of (model, Z, losses) - only valid on rank 0
        Other ranks return (None, None, [])

    Example:
        >>> # Run with torchrun
        >>> # torchrun --nproc_per_node=128 script.py
        >>> model, Z, losses = train_autoencoder_ddp(X, config)

    Usage (ARCHER2):
        ```bash
        torchrun --nproc_per_node=128 \\
            bin/train_ddp.py --aa L --encoding aa_20 --distance 15
        ```
    """
    # Initialize DDP
    rank = _get_rank()
    world_size = _get_world_size()

    if world_size > 1:
        # Use gloo backend for CPU, nccl for GPU
        backend = "nccl" if torch.cuda.is_available() else "gloo"
        dist.init_process_group(backend=backend)

    # Set device
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{rank % torch.cuda.device_count()}")
    else:
        device = torch.device("cpu")

    # Only rank 0 logs
    verbose = (rank == 0)

    if verbose:
        logger.info(
            f"Training {model_type} autoencoder (DDP): "
            f"{X.shape[1]}D → {config.latent_dim}D "
            f"({X.shape[0]:,} samples, {config.epochs} epochs, {world_size} processes)"
        )

    # Validate input
    if X.ndim != 2:
        raise ValueError(f"X must be 2D, got shape {X.shape}")
    if X.dtype != np.float32:
        X = X.astype(np.float32)

    n_samples, input_dim = X.shape

    # Create model
    if model_type == "standard":
        model = Autoencoder(input_dim=input_dim, latent_dim=config.latent_dim)
    elif model_type == "variational":
        model = VariationalAutoencoder(
            input_dim=input_dim,
            latent_dim=config.latent_dim,
            beta=0.1,
        )
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    model = model.to(device)

    if verbose:
        logger.info(f"  Device: {device}")
        logger.info(f"  Model: {model.count_parameters():,} parameters")
        logger.info(f"  World size: {world_size}")

    # Wrap with DDP
    if world_size > 1:
        model = nn.parallel.DistributedDataParallel(
            model,
            device_ids=[rank] if torch.cuda.is_available() else None,
            find_unused_parameters=False,
        )

    # Distributed data sampler (critical for DDP!)
    X_tensor = torch.FloatTensor(X).to(device)
    dataset = TensorDataset(X_tensor)

    sampler = DistributedSampler(
        dataset,
        num_replicas=world_size,
        rank=rank,
        shuffle=True,
        seed=42,  # For reproducibility
    )

    loader = DataLoader(
        dataset,
        batch_size=config.batch_size,
        sampler=sampler,
        pin_memory=(device.type == "cuda"),
    )

    # Optimizer
    if config.optimizer.lower() == "adam":
        optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)
    elif config.optimizer.lower() == "sgd":
        optimizer = torch.optim.SGD(model.parameters(), lr=config.learning_rate)
    else:
        raise ValueError(f"Unknown optimizer: {config.optimizer}")

    # Loss function
    if config.loss_function.lower() == "mse":
        loss_fn = nn.MSELoss(reduction="mean")
    elif config.loss_function.lower() == "bce":
        loss_fn = nn.BCELoss(reduction="mean")
    else:
        raise ValueError(f"Unknown loss: {config.loss_function}")

    if verbose:
        logger.info(f"  Optimizer: {config.optimizer}, lr={config.learning_rate}")
        logger.info(f"  Loss: {config.loss_function}")

    # Training loop
    losses = []
    for epoch in range(config.epochs):
        sampler.set_epoch(epoch)  # Important: reshuffle per epoch in DDP
        model.train()
        epoch_losses = []

        for (batch_X,) in loader:
            optimizer.zero_grad()

            if model_type == "standard":
                batch_recon = model(batch_X)
                loss = loss_fn(batch_recon, batch_X)
            else:
                batch_recon, mean, log_var = model(batch_X)
                recon_loss = loss_fn(batch_recon, batch_X)
                kl_loss = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
                kl_loss = kl_loss / batch_X.shape[0]
                loss = recon_loss + model.beta * kl_loss

            loss.backward()
            optimizer.step()
            epoch_losses.append(loss.item())

        epoch_loss = np.mean(epoch_losses)
        losses.append(epoch_loss)

        # Periodic logging (rank 0 only)
        if verbose and (epoch + 1) % max(1, config.epochs // 5) == 0:
            logger.info(f"  Epoch {epoch + 1:3d}/{config.epochs}: loss={epoch_loss:.4f}")

    if verbose:
        logger.info(f"✓ Training complete: final loss={losses[-1]:.4f}")

    # Extract embeddings (rank 0 only)
    if rank == 0:
        model_unwrapped = model.module if isinstance(model, nn.parallel.DistributedDataParallel) else model
        model_unwrapped.eval()

        with torch.no_grad():
            if model_type == "standard":
                Z_tensor = model_unwrapped.encode(X_tensor)
            else:
                mean, _ = model_unwrapped.encode(X_tensor)
                Z_tensor = mean

        Z = Z_tensor.cpu().numpy().astype(np.float32)
        model_unwrapped = model_unwrapped.cpu()
        return model_unwrapped, Z, losses
    else:
        # Non-rank-0 processes return empty results
        return None, np.array([]), []


def train_autoencoder_with_validation(
    X_train: np.ndarray,
    X_val: np.ndarray,
    config: TrainingConfig,
    model_type: str = "standard",
    verbose: bool = True,
) -> Tuple[nn.Module, np.ndarray, np.ndarray, Dict]:
    """Train autoencoder with validation set.

    Wrapper around train_autoencoder that handles train/val split and returns
    both training and validation embeddings.

    Args:
        X_train: Training data (n_train, input_dim)
        X_val: Validation data (n_val, input_dim)
        config: TrainingConfig instance
        model_type: Type of autoencoder
        verbose: Whether to log progress

    Returns:
        Tuple of (model, Z_train, Z_val, history) where history contains:
            - 'train_losses': Per-epoch training losses
            - 'val_losses': Per-epoch validation losses (if early stopping enabled)

    Example:
        >>> X_train = np.random.randn(8000, 900).astype(np.float32)
        >>> X_val = np.random.randn(2000, 900).astype(np.float32)
        >>> config = TrainingConfig(latent_dim=16, epochs=10, early_stopping_patience=3)
        >>> model, Z_train, Z_val, hist = train_autoencoder_with_validation(X_train, X_val, config)
    """
    model, Z_train, losses = train_autoencoder(
        X_train,
        config,
        model_type=model_type,
        validate_x=X_val,
        verbose=verbose,
    )

    # Extract validation embeddings
    X_val_tensor = torch.FloatTensor(X_val).to(torch.device(config.device))
    model = model.to(torch.device(config.device))
    model.eval()

    with torch.no_grad():
        if model_type == "standard":
            Z_val_tensor = model.encode(X_val_tensor)
        else:
            mean, _ = model.encode(X_val_tensor)
            Z_val_tensor = mean

    Z_val = Z_val_tensor.cpu().numpy().astype(np.float32)
    model = model.cpu()

    history = {
        "train_losses": losses,
    }

    return model, Z_train, Z_val, history


def save_model(model: nn.Module, path: str, config: Optional[dict] = None) -> None:
    """Save trained autoencoder to disk.

    Args:
        model: Trained model instance
        path: Path to save .pt file
        config: Optional config dict to save alongside model

    Example:
        >>> save_model(model, "autoencoder.pt", config=model.get_config())
    """
    torch.save(model.state_dict(), path)
    logger.info(f"Saved model: {path}")

    if config:
        import json
        config_path = path.replace(".pt", "_config.json")
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)
        logger.info(f"Saved config: {config_path}")


def load_model(model: nn.Module, path: str) -> nn.Module:
    """Load trained autoencoder from disk.

    Args:
        model: Initialized model instance
        path: Path to .pt file

    Returns:
        Model with loaded weights

    Example:
        >>> model = Autoencoder(900, 16)
        >>> model = load_model(model, "autoencoder.pt")
    """
    model.load_state_dict(torch.load(path))
    logger.info(f"Loaded model: {path}")
    return model
