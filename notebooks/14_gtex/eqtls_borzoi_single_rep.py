import os
import grelu.resources
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
from decima.data.dataset import VariantDataset


def undo_squashed_scale(x, clip_soft=384, track_transform=3 / 4):
    x = x.clone()
    unclip_mask = x > clip_soft
    if unclip_mask.any():
        x[unclip_mask] = (x[unclip_mask] - clip_soft) ** 2 + clip_soft

    x = (x + 1) ** (1.0 / track_transform) - 1
    return x


parser = argparse.ArgumentParser()
parser.add_argument("--variants", type=str)
parser.add_argument("--out_dir", type=str)
parser.add_argument("--prefix", type=str)
parser.add_argument("--device", type=str)
parser.add_argument("--rep", type=int)
parser.add_argument("--chunkstart", type=int, default=None)
parser.add_argument("--chunkend", type=int, default=None)
args = parser.parse_args()


os.environ["CUDA_VISIBLE_DEVICES"] = args.device
import torch
from torch.utils.data import DataLoader

print("Loading variants")
variants = pd.read_csv(args.variants)
n_variants = len(variants)
print(f"Read {n_variants} variants")

if args.chunkstart is not None:
    print(f"Subsetting variants: {args.chunkstart}-{args.chunkend}")
    variants = variants.iloc[args.chunkstart:args.chunkend]
    n_variants = len(variants)
    print(f"Subsetted to {n_variants} variants")


print("Creating dataset")
ds = VariantDataset(
    variants,
    metadata_anndata="/gstore/data/resbioai/grelu/decima/20240823/data.h5ad",
    gene_col="gene",
)
assert len(ds.variants)==n_variants


print("Making bins")
bins = pd.DataFrame(
    {
        "start": ds.variants.gene_mask_start.apply(lambda x: int(np.floor(((x) / 32) - 5120))),
        "end": ds.variants.gene_mask_end.apply(lambda x: int(np.floor(((x) / 32) - 5120))) + 1,
    }
)

print("Generating predictions")
device = torch.device(0)

with torch.no_grad():
    print(f"Replicate: {args.rep}")
    print("Loading model")
    model = grelu.resources.load_model(
        project="borzoi",
        model_name=f"human_rep{args.rep}",
    )
    model.model.eval()
    model = model.to(device)

    print("Creating dataloader")
    dl = DataLoader(ds, batch_size=2, shuffle=False)
    dl = iter(dl)

    print("Making predictions")
    preds = []
    for i in tqdm(range(len(variants))):
        bin = bins.iloc[i]
        seqs = next(dl)['seq'][:, :4]
        seqs = seqs.to(device)
        
        pred = model(seqs).detach().cpu()  # 2, 7611, 6144
        pred = undo_squashed_scale(pred)

        # Sum predicted coverage over gene
        pred = pred[:, :, bin.start: bin.end + 1]  # 2, 7611, l
        pred = pred.sum(-1).numpy()  # 2, 7611

        # Log total coverage
        pred = np.log(pred + 1)  # 2, 7611

        # Alt - ref
        pred = pred[1] - pred[0]  # 7611
        preds.append(pred)


print("Combining predictions")
preds = np.stack(preds) # len(variants), 7611
print(f"Predictions shape: {preds.shape}") # len(variants), 7611

print("Saving predictions")
if args.chunkstart is None:
    out_file = os.path.join(args.out_dir, f"{args.prefix}_borzoi_{args.rep}.npy")
else:
    out_file = os.path.join(args.out_dir, f"{args.prefix}_borzoi_{args.rep}_{args.chunkstart}_{args.chunkend}.npy")
np.save(out_file, preds)
