#!/bin/bash
#SBATCH --account=owner_gred_braid_gpu
#SBATCH --partition=owner_gred_braid_gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --job-name=decima_gtex       # Job name for the array
#SBATCH --output=decima_gtex_%A_%a.out  # %A for array job ID, %a for array task ID
#SBATCH --error=decima_gtex_%A_%a.err   #
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --array=0-3

# --- SLURM Job Array Configuration ---
reps=(0 1 2 3)

# Load necessary modules 
module load CEDAR/2020.08
module load Python/gpy_gpu/Python310
python3 -m pip install git+https://github.com/Genentech/decima.git
python3 -m pip install --upgrade grelu
export PATH="/home/lala8/.local/bin:$PATH"

script="/home/lala8/decima-applications/avantika_notebooks/20240823/14_bulk_eqtls/eqtls_decima_single_rep.py"
variants="/home/lala8/decima-applications/avantika_notebooks/20240823/14_bulk_eqtls/gtex_eqtl_cat/variants_df.csv"

python3 $script --variants $variants --out_dir /gstore/data/resbioai/grelu/decima/20240823/bulk_eqtl_results --prefix gtex_eqtl_cat --device 0 --rep $SLURM_ARRAY_TASK_ID