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
parser.add_argument("-fold", 
                    help="Which brozoi fold to use",
                    type=int)
parser.add_argument("-gene_h5_file", 
                    help="Path to h5 file")
parser.add_argument("-variant_df_file", 
                    help="Path to variant dataframe")
parser.add_argument("-out_dir", 
                    help="Directory where to save results")
parser.add_argument("-unsquash", 
                    help="Whether to unsquash", action='store_true')
parser.add_argument("-tracks", 
                    help="Path to borzoi tracks")
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
print(start_idx, end_idx)

fold = args.fold
unsquash = args.unsquash
tracks_file = args.tracks
print(unsquash)

gene_h5_file = args.gene_h5_file
variant_df_file = args.variant_df_file
out_dir=args.out_dir
suffix=args.suffix
unsquash_suffix = "_unsquash" if unsquash else ""
if not os.path.exists(os.path.join(out_dir,f'brozima_fold{fold}'+suffix+unsquash_suffix)):
    os.mkdir(os.path.join(out_dir,f'brozima_fold{fold}'+suffix+unsquash_suffix))
out_dir_gene=os.path.join(out_dir,f'brozima_fold{fold}'+suffix+unsquash_suffix,f'gene_scores_{start_idx}_{end_idx}.npy')
out_dir_tss=os.path.join(out_dir,f'brozima_fold{fold}'+suffix+unsquash_suffix,f'tss_scores_{start_idx}_{end_idx}.npy')
print(out_dir_gene)

seq_len = 524288
batch_size = 12

variant_df = pd.read_csv(variant_df_file)
variant_df = variant_df.iloc[start_idx:end_idx]
print(len(variant_df))

# Import here after setting visible devices
from grelu.sequence.format import strings_to_indices, indices_to_one_hot
from grelu.sequence.mutate import mutate
from grelu.sequence.utils import reverse_complement

import torch
from torch import nn

from torch.utils.data import DataLoader, Dataset

# Dataset for variant effect prediction
class eQTLDataset(Dataset):
    def __init__(self, seq_h5_file, eqtl_df, 
                 return_ref_check=False, rc=False, seq_len=524288):
        self.h5_file = seq_h5_file

        with h5py.File(self.h5_file, 'r') as f:
            genes = np.array(f["genes"]).astype(str)
            genes = genes[:, 0]
        self.gene_ord = genes
   
        self.eqtl_df = eqtl_df # this must have gene_symbol, variant, chrom, ref, alt, pos_relative
        self.eqtl_len = len(self.eqtl_df)
        self.dataset = None 
        self.return_ref_check = return_ref_check

        self.n_augmented = 1
        self.n_alleles = 1
        self.rc = rc
        self.seq_len = seq_len

    def __len__(self):
        return self.eqtl_len * (4 if self.rc else 2)


    def _extract_center(self, x, shift=0):
        start = (x.shape[-1] - self.seq_len)//2
        start -= shift
        return x[..., start:start+self.seq_len]

    def extract_seq(self, idx):
        seq = self.dataset['sequences'][idx]
        seq = indices_to_one_hot(seq) # 4, L
        mask = self.dataset['masks'][[idx]] # 1, L
        # extract center
        seq = self._extract_center(seq)
        mask = self._extract_center(mask)
        # combine
        seq = np.concatenate([seq, mask]) # 5, L
        return torch.Tensor(seq)

    def __getitem__(self, idx):
        if idx >= self.__len__():
            raise IndexError

        if self.dataset is None:
            self.dataset = h5py.File(self.h5_file, 'r')
        
        # select rc based on index, second half is rc e.g. ref/alt/ref-rc/alt-rc
        if self.rc and idx >= (self.eqtl_len*2):
            rc = True
            idx = idx % (self.eqtl_len*2)
        else:
            rc = False

        # split into ref sequences first, alt sequences next
        allele = 'ref' if idx < self.eqtl_len else 'alt'
        rec = self.eqtl_df.iloc[idx % self.eqtl_len]

        # get eQTL info
        gene = rec.gene_symbol
        ref = strings_to_indices(rec.ref)[0]
        alt = strings_to_indices(rec.alt)[0]
        offset = rec.pos_relative
        
        # Get sequence and mask
        seq_idx = np.where(self.gene_ord == gene)[0][0]
        inp = self.extract_seq(seq_idx)
        #print(inp.shape)
        seq = inp[:4].argmax(axis=0) # convert to indices
        mask = inp[4]

        assert seq[offset] == ref, rec.variant + "_" + rec.gene_symbol

        seq = mutate(seq, allele=ref if allele == 'ref' else alt, pos=offset, input_type='indices')
        if rc:
            seq = reverse_complement(seq, input_type='indices')

        seq = indices_to_one_hot(seq)
        seq = torch.Tensor(seq)

        # pool mask
        mask = torch.nn.functional.max_pool1d(mask.unsqueeze(0),32).squeeze()

        # compute tss mask
        tss_pos = mask.argmax()
        tss_mask = torch.zeros_like(mask)
        tss_mask[tss_pos - 10: tss_pos + 10] = 1

        mask = mask.float().unsqueeze(0)
        tss_mask = tss_mask.float().unsqueeze(0)
        return seq,mask,tss_mask

# Load model

import grelu.resources
model = grelu.resources.load_model(
    project="borzoi",
    model_name=f"human_fold{fold}",
)
model = model.model
# remove cropping (as borzoi paper does it too for eQTL)
model.embedding.crop = torch.nn.Identity()
model.eval()
model.to(device)

# Construct unsquashing tensors
# (Since sc-data is in log1p, not squashed scale)
brozoi_tracks = pd.read_csv(tracks_file,sep='\t')
clip_soft = torch.tensor(np.array(brozoi_tracks['clip_soft']))
clip_soft = clip_soft.unsqueeze(-1).expand(-1,16384).unsqueeze(0).to(device)
track_scale = torch.tensor(np.array(brozoi_tracks['scale']))
track_scale = track_scale.unsqueeze(-1).expand(-1,16384).unsqueeze(0).to(device)
track_transform = 3/4

# Run
var_ds = eQTLDataset(gene_h5_file, variant_df, rc=False)
var_dl = DataLoader(var_ds, shuffle=False, batch_size=12, num_workers=8, pin_memory=True)

gene_preds = []
tss_preds = []

with torch.no_grad():
    for batch in tqdm.tqdm(var_dl):
        seq, mask, tss_mask = batch
        seq, mask, tss_mask = seq.to(device), mask.to(device), tss_mask.to(device)
        pred = model(seq)
        # gene-level prediction:
        # select gene bins with a gene-mask
        # Sum across the gene mask
        if unsquash:
            # undo scale
            gene_pred = pred / track_scale
            # undo soft_clip
            if clip_soft is not None:
                gene_pred = torch.where(
                    gene_pred > clip_soft, (gene_pred - clip_soft) ** 2 + clip_soft, gene_pred
                )
            # undo sqrt
            gene_pred = gene_pred ** (1.0 / track_transform)
            # mask and sum
            gene_pred = gene_pred*mask
            gene_pred = gene_pred.sum(axis=-1)
        else:
            masked_pred = pred*mask
            gene_pred = masked_pred.sum(axis=-1)
        # TSS-level prediction (for CAGE):
        tss_pred = (pred*tss_mask).sum(axis=-1)
        # move to cpu
        gene_pred = gene_pred.cpu()
        tss_pred = tss_pred.cpu()
        gene_preds.append(gene_pred)
        tss_preds.append(tss_pred)

# Save results

gene_preds = torch.concat(gene_preds,axis=0).numpy()
tss_preds = torch.concat(tss_preds,axis=0).numpy()

np.save(out_dir_gene, gene_preds)
np.save(out_dir_tss, tss_preds)





