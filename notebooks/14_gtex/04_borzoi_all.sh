#!/bin/bash
#SBATCH --account=owner_gred_braid_gpu
#SBATCH --partition=owner_gred_braid_gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --job-name=borzoi_gtex       # Job name for the array
#SBATCH --output=borzoi_gtex_%A_%a.out  # %A for array job ID, %a for array task ID
#SBATCH --error=borzoi_gtex_%A_%a.err   #
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --array=0-47

# --- SLURM Job Array Configuration ---
reps=(0 1 2 3)
starts=(0 20000 40000 60000 80000 100000 120000 140000 160000 180000 200000 220000)
ends=(20000 40000 60000 80000 100000 120000 140000 160000 180000 200000 220000 229828)

# Calculate total combinations for the array size
num_reps=${#reps[@]}
num_chunks=${#starts[@]}

# --- Hyperparameter Selection for the Current Array Task ---
rep_index=$((SLURM_ARRAY_TASK_ID % num_reps))
chunk_index=$((SLURM_ARRAY_TASK_ID / num_reps % num_chunks))

REP=${reps[$rep_index]}
START=${starts[$chunk_index]}
END=${ends[$chunk_index]}

# Load necessary modules 
module load CEDAR/2020.08
module load Python/gpy_gpu/Python310
python3 -m pip install git+https://github.com/Genentech/decima.git
python3 -m pip install --upgrade grelu
export PATH="/home/lala8/.local/bin:$PATH"

script="/home/lala8/decima-applications/avantika_notebooks/20240823/14_bulk_eqtls/eqtls_borzoi_single_rep.py"
variants="/gstore/data/resbioai/grelu/decima/20240823/bulk_eqtl_results/variants_df.csv"

python3 $script --variants $variants --out_dir /gstore/data/resbioai/grelu/decima/20240823/bulk_eqtl_results --prefix gtex_eqtl_cat --device 0 --rep $REP --chunkstart $START --chunkend $END
