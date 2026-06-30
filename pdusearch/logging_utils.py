"""Centralized logging utilities for PDU extraction and analysis."""

import logging
import sys
from pathlib import Path
from typing import Optional
from datetime import datetime


class ColoredFormatter(logging.Formatter):
    """Colored formatter for console output."""

    COLORS = {
        "DEBUG": "\033[36m",      # Cyan
        "INFO": "\033[32m",       # Green
        "WARNING": "\033[33m",    # Yellow
        "ERROR": "\033[31m",      # Red
        "CRITICAL": "\033[35m",   # Magenta
    }
    RESET = "\033[0m"

    def format(self, record):
        """Format log record with colors."""
        levelname = record.levelname
        color = self.COLORS.get(levelname, self.RESET)
        record.levelname = f"{color}{levelname}{self.RESET}"
        return super().format(record)


def configure_logging(
    log_file: Optional[str] = None,
    log_level: str = "INFO",
    name: str = "pdusearch",
    use_colors: bool = True,
) -> logging.Logger:
    """Configure logging for PDU extraction and analysis.

    Args:
        log_file: Optional path to log file. If None, only console logging.
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        name: Logger name.
        use_colors: Whether to use colored console output.

    Returns:
        Configured logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, log_level.upper()))

    # Remove existing handlers to avoid duplication
    logger.handlers.clear()

    # Create formatters
    if use_colors:
        console_formatter = ColoredFormatter(
            fmt="%(levelname)-8s [%(name)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    else:
        console_formatter = logging.Formatter(
            fmt="%(levelname)-8s [%(name)s] %(asctime)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    file_formatter = logging.Formatter(
        fmt="%(levelname)-8s [%(name)s] %(asctime)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, log_level.upper()))
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setLevel(getattr(logging, log_level.upper()))
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger


def add_logging_args(parser):
    """Add standard logging arguments to ArgumentParser.

    Args:
        parser: argparse.ArgumentParser instance.
    """
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--log-file",
        default=None,
        help="Log file path (default: None, logs to console only)",
    )


def get_logger(name: str = "pdusearch") -> logging.Logger:
    """Get or create a logger.

    Args:
        name: Logger name.

    Returns:
        Logger instance.
    """
    return logging.getLogger(name)


def log_timestamp(logger: logging.Logger, message: str):
    """Log a message with timestamp.

    Args:
        logger: Logger instance.
        message: Message to log.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logger.info(f"{message} [{timestamp}]")


def log_separator(logger: logging.Logger, char: str = "=", width: int = 80):
    """Log a separator line.

    Args:
        logger: Logger instance.
        char: Character to use for separator.
        width: Width of separator.
    """
    logger.info(char * width)
