import subprocess
import sys
import math
import numpy as np
import pandas as pd
import anndata
import h5py
import os
import scipy
import tqdm
import scanpy as sc
import wandb

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-device", 
                    help="which gpu to use",
                    type=int)
parser.add_argument("-task", 
                    help="Which chunk should we predict",
                    type=int)
parser.add_argument("-job_size", 
                    help="How big is each chunk",
                    type=int)
parser.add_argument("-ckpt_file", 
                    help="Path to the model checkpoint")
parser.add_argument("-gene_h5_file", 
                    help="Path to h5 file")
parser.add_argument("-variant_df_file", 
                    help="Path to variant dataframe")
parser.add_argument("-out_dir", 
                    help="Directory where to save results")
parser.add_argument("-suffix", 
                    help="Suffix to output directory",
                    default='')
args = parser.parse_args()
device = args.device
os.environ["CUDA_VISIBLE_DEVICES"] = str(device)
device = 0

task = args.task
job_size = args.job_size
start_idx = (task * job_size)
end_idx = ((task+1) * job_size)

ckpt_file=args.ckpt_file

suffix = args.suffix
gene_h5_file = args.gene_h5_file
variant_df_file = args.variant_df_file
out_dir=args.out_dir
if not os.path.exists(os.path.join(out_dir,'eqtl'+suffix)):
    os.mkdir(os.path.join(out_dir,'eqtl'+suffix))
out_dir=os.path.join(out_dir,'eqtl'+suffix,f'eqtl_scores_{start_idx}_{end_idx}.npy')

batch_size = 12

variant_df = pd.read_csv(variant_df_file)
variant_df = variant_df.iloc[start_idx:end_idx]

sys.path.insert(1, '/home/karollua/projects/Decima/scborzoi/decima/src')
from decima.read_hdf5 import VariantDataset
from decima.lightning import LightningModel

# Load model
model = LightningModel.load_from_checkpoint(ckpt_file)
model = model.eval()

# Run
var_ds = VariantDataset(variant_df.rename(columns={'gene_symbol':'gene','pos_relative':'rel_pos','ref':'ref_tx','alt':'alt_tx'}),gene_h5_file,test_ref=True)
scores = model.predict_on_dataset(var_ds, devices=device, batch_size=12, num_workers=8).squeeze()

# Save results
np.save(out_dir, scores)

print("All done")




