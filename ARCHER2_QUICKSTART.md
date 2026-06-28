# ARCHER2 Quick Start Guide

Account: `e280-dantu`

## Step 1: Login to ARCHER2

```bash
ssh username@login.archer2.ac.uk
```

## Step 2: Transfer Code and Features

From your local machine:

```bash
# Copy pdu-extractor code
scp -r /Users/sarathdantu/Projects/vsprojects/mdpdu/pdu-extractor \
    username@login.archer2.ac.uk:~/pdu-extractor

# Copy features (or just features for the AAs you need)
scp -r /path/to/analysis/features \
    username@login.archer2.ac.uk:~/pdu-extractor/analysis/
```

Or if you're already on ARCHER2 and have shared storage access, copy from there.

## Step 3: Setup Virtual Environment (One-Time)

On ARCHER2 login node:

```bash
cd ~
bash pdu-extractor/setup_archer2_venv.sh
```

This creates `venv_archer2/` with all dependencies.

**Output should show:**
```
✓ Virtual environment ready
```

## Step 4: Submit Training Jobs

### Single AA (LEU), Single Node

```bash
sbatch pdu-extractor/submit_archer2.slurm
```

Expected runtime: **~2-3 minutes** for 2.1M PDUs

### Single AA (different amino acid)

```bash
sbatch -J PDU_TRAIN_A --export=PDU_AA=A pdu-extractor/submit_archer2.slurm
```

### Multiple Nodes (4 nodes = 512 cores)

```bash
sbatch -N 4 pdu-extractor/submit_archer2.slurm
```

Expected runtime: **~30-45 seconds**

### All 20 AAs in Parallel

Create `train_all.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=PDU_ALL_AA
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --array=0-19

module load cray-python
source venv_archer2/bin/activate

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

AAS=(A C D E F G H I K L M N P Q R S T V W Y)
AA=${AAS[$SLURM_ARRAY_TASK_ID]}

mkdir -p logs

srun python3 pdu-extractor/scripts/train_pdu_autoencoder_archer2.py \
    --features-dir pdu-extractor/analysis/features \
    --out-dir pdu-extractor/analysis/embeddings \
    --aa "$AA" \
    --batch-size 1024 \
    --epochs 15
```

Submit:
```bash
sbatch train_all.slurm
```

This submits 20 jobs (one per AA) running in parallel on 2 nodes.

## Step 5: Monitor Jobs

```bash
# List your jobs
squeue -u $USER

# Watch specific job
watch squeue -j <jobid>

# Check completion
squeue -j <jobid> --format=JobID,JobName,State,Elapsed

# Check if finished
sacct -j <jobid> --format=JobID,JobName,State,ExitCode
```

## Step 6: Check Output

While job is running:

```bash
tail -f logs/pdu_train_<jobid>.out
```

Expected output:
```
ARCHER2 PDU Autoencoder Training
AA: L
Nodes: 1
Tasks: 128
...

Training on ARCHER2: 1 processes, rank 0
Loaded 2,086,440 PDU features (shape: (2086440, 900))

Epoch 1/15, Batch 100/2041, Loss: 0.000234
Epoch 1/15, Batch 200/2041, Loss: 0.000198
...

Epoch 15/15: train_loss=0.000184, val_loss=0.000182
✓ Training complete
```

## Step 7: Copy Results Back

After training completes, copy embeddings to your Mac:

```bash
scp -r username@login.archer2.ac.uk:~/pdu-extractor/analysis/embeddings \
    /Users/sarathdantu/Projects/vsprojects/mdpdu/pdu-extractor/analysis/
```

## Batch Job Parameters

In `submit_archer2.slurm`:

- `--nodes=1` → change to `-N 4` for 4 nodes
- `--ntasks=128` → total tasks (128 = 1 node)
- `--tasks-per-node=128` → tasks per node (always 128)
- `--cpus-per-task=1` → keep at 1 (full core per task)
- `--time=2:00:00` → increase for more complex training

Memory per node: 256 GB (standard) → batch_size 1024 is safe

## Troubleshooting

### Job fails: "Virtual environment not found"

```bash
bash pdu-extractor/setup_archer2_venv.sh
```

### Job fails: "No such file or directory"

Check paths are correct:
```bash
ls -la pdu-extractor/analysis/features/pdu_features_L.npz
```

### Job slow or not scaling

Check node allocation:
```bash
sinfo -n <nodename> -l
```

Check CPU affinity:
```bash
srun -n 1 numactl -s
```

### Out of memory

Reduce batch size in submit script:
```bash
--batch-size 512
```

## Cost Estimate

ARCHER2 charges in "core-hours":
- 1 node × 2 hours = 128 × 2 = 256 core-hours
- 4 nodes × 0.5 hours = 512 × 0.5 = 256 core-hours

Rough cost: £1-2 per job

## Next Steps

1. All embeddings trained ✓
2. Copy to local Mac ✓
3. Run HDBSCAN clustering locally
4. Run PyCoM enrichment analysis
5. Check if 740 LEU clusters have biological meaning

## Support

ARCHER2 support: support@archer2.ac.uk
Account issues: account@archer2.ac.uk
