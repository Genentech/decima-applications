import numpy as np
import pandas as pd
import anndata
import h5py
import os
import scipy
import tqdm
import gc
import sys

# Parameters

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-device", 
                    help="Which gpu to use",
                    type=int)
parser.add_argument("-ckpt_file", 
                    help="Path to the model checkpoint")
parser.add_argument("-h5_file", 
                    help="Path to h5 file (gene-indexed)")
parser.add_argument("-gene_df_file", 
                    help="Path to gene df (for which to attribute)")
parser.add_argument("-targets_file", 
                    help="Path to df containing on and off tasks")
parser.add_argument("-out_dir", 
                    help="Directory where to save results")
parser.add_argument("-transform", 
                    help="Which kind of attributions to compute",
                    default="specificity")
args = parser.parse_args()
device = args.device
os.environ["CUDA_VISIBLE_DEVICES"] = str(device)
device = 0

gene_df_file = args.gene_df_file
targets_file = args.targets_file
h5_file = args.h5_file
ckpt_file = args.ckpt_file
out_dir = args.out_dir
transform = args.transform

# ipmort here after setting visible devices
from captum.attr import Saliency
sys.path.insert(1, '/home/karollua/projects/Decima/scborzoi/decima/src')
from decima.lightning import LightningModel
from decima.interpret import attributions as get_attr

# Load and Prepare data
gene_df = pd.read_csv(gene_df_file)
targets_df = pd.read_csv(targets_file)

# Get tasks and genes
target_genes = gene_df['gene'].tolist()
on_tasks = targets_df.query('task_type == "on"')['task'].tolist()
off_tasks = targets_df.query('task_type == "off"')['task'].tolist()
if len(off_tasks) == 0:
    off_tasks = None

# Prepare model
model = LightningModel.load_from_checkpoint(ckpt_file)
model = model.eval()

# Run Attribution
print('Attributing')
sequences = []
attributions = []
for gene in tqdm.tqdm(target_genes):
    seq, tss_pos, attr = get_attr(gene=gene, 
                                  h5_file=h5_file, 
                                  model=model, 
                                  device=device, 
                                  tasks=on_tasks,
                                  off_tasks=off_tasks,
                                  transform=transform,
                                  method=Saliency,
                                  abs=False)

    sequences.append(seq)
    attributions.append(attr)

attributions = np.stack(attributions)
sequences = np.stack(sequences)

# Save attributions
np.save(os.path.join(out_dir, 'attributions.npy'),attributions)
np.save(os.path.join(out_dir, 'sequences.npy'),sequences)

print("All done")