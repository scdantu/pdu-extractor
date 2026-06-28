# ARCHER2 Distributed Autoencoder Training

## Overview

Train PDU autoencoders on ARCHER2's CPU nodes using distributed training (PyTorch DDP).

**Performance Estimates:**
- Local Mac (8 cores): 517K PDUs in 47 min = ~180 PDUs/sec
- ARCHER2 1 node (128 cores): **~20-30x faster** = ~3600-5400 PDUs/sec
- ARCHER2 4 nodes (512 cores): **~60-100x faster** = ~10K-18K PDUs/sec

**For all 20 AAs:**
- Local Mac: ~30 hours
- ARCHER2 (1 node): ~1-2 hours
- ARCHER2 (4 nodes): ~20-30 minutes

## Setup

### 1. Connect to ARCHER2

```bash
ssh username@login.archer2.ac.uk
```

### 2. Transfer Code & Data

```bash
# Copy pdu-extractor to ARCHER2
scp -r /Users/sarathdantu/Projects/vsprojects/mdpdu/pdu-extractor \
    username@login.archer2.ac.uk:~/pdu-extractor

# Copy features (if not already on ARCHER2)
# Features are large; consider copying only what you need per AA
scp -r /path/to/features username@login.archer2.ac.uk:~/pdu-extractor/analysis/features
```

### 3. Check Account

Edit `submit_archer2.slurm` and set your ARCHER2 account:

```bash
#SBATCH --account=<YOUR_ACCOUNT>
```

Get your account:
```bash
sacctmgr show user $USER --format=User,DefaultAccount
```

## Training: Single AA

### Option 1: Single Node (128 cores)

```bash
sbatch -J pdu_L pdu-extractor/submit_archer2.slurm L
```

**Expected runtime:** ~2-3 minutes for 2.1M PDUs (LEU)

### Option 2: Multiple Nodes (4 nodes = 512 cores)

```bash
sbatch -N 4 -J pdu_L pdu-extractor/submit_archer2.slurm L
```

**Expected runtime:** ~30-45 seconds for 2.1M PDUs

## Training: All 20 AAs (Batch Job)

Create `train_all_archer2.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=pdu_all_aa
#SBATCH --nodes=2
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --array=0-19

module load PrgEnv-gnu
module load cray-python

# Map array index to amino acid
AAS=(A C D E F G H I K L M N P Q R S T V W Y)
AA=${AAS[$SLURM_ARRAY_TASK_ID]}

echo "Training AA=$AA (task $SLURM_ARRAY_TASK_ID)"

python3 pdu-extractor/scripts/train_pdu_autoencoder_archer2.py \
    --features-dir pdu-extractor/analysis/features \
    --out-dir pdu-extractor/analysis/embeddings \
    --aa "$AA" \
    --batch-size 1024 \
    --epochs 15
```

Submit:
```bash
sbatch train_all_archer2.slurm
```

This submits 20 parallel jobs (one per AA) on 2 nodes.

## Monitor Jobs

```bash
# List your jobs
squeue -u $USER

# Check specific job
squeue -j <jobid>

# Check job details
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,CPUTime,MaxRSS

# Watch real-time output
tail -f logs/pdu_train_<jobid>.out
```

## Performance Tuning

### Batch Size

In `submit_archer2.slurm`, adjust `--batch-size`:
- Single node: 512-1024 (256GB RAM)
- Multi-node: 2048-4096 (distributed across nodes)

Larger batches = faster training but more memory.

### Threads

In submit script, adjust MKL threads:
```bash
export MKL_NUM_THREADS=4  # Threads per PyTorch worker
export OMP_NUM_THREADS=4
```

More threads = better parallelization but can oversub scribe. Rule: `threads × workers ≈ CPUs / 2`

### CPU Affinity

```bash
export KMP_AFFINITY=balanced  # NUMA-aware affinity
```

## Expected Output

```
Training on ARCHER2: 2 processes, rank 0
Batch size: 1024, Epochs: 15
Loaded 2,086,440 PDU features (shape: (2086440, 900))

Epoch 1/15, Batch 100/2041, Loss: 0.000234
Epoch 1/15, Batch 200/2041, Loss: 0.000198
...
Epoch 15/15: train_loss=0.000184, val_loss=0.000182

Saved embeddings: analysis/embeddings/pdu_embedding_L.npz
Shape: (2086440, 16)
Best validation loss: 0.000182
```

## Troubleshooting

### Job fails with "Module not found"

Check available Python modules:
```bash
module avail cray-python
module load cray-python  # Load latest
python3 -m pip list | grep torch
```

### Out of memory

Reduce batch size:
```bash
--batch-size 512  # Instead of 1024
```

### Slow performance

Check CPU affinity:
```bash
srun -n 1 -c 128 numactl -s
```

Contact ARCHER2 support: support@archer2.ac.uk

## Copy Results Back

After training, copy embeddings and clusters back to your Mac:

```bash
scp -r username@login.archer2.ac.uk:~/pdu-extractor/analysis/embeddings \
    /Users/sarathdantu/Projects/vsprojects/mdpdu/pdu-extractor/analysis/
```

## Cost Estimation

ARCHER2 charges in Core-Hours (node × hours):
- 1 node × 2 hours = 128 core-hours ≈ £0.50 (rough estimate)
- 20 AAs on 2 nodes × 1 hour = ~£50-100 total budget-friendly

## Next Steps

After embeddings are ready on ARCHER2:

1. **Copy embeddings to local Mac**
2. **Run clustering locally** (faster HDBSCAN on recent embeddings)
3. **Run PyCoM enrichment** to analyze if 740 clusters have biological meaning
