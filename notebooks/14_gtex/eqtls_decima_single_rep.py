import os
import anndata
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--variants", type=str)
parser.add_argument("--out_dir", type=str)
parser.add_argument("--prefix", type=str)
parser.add_argument("--device", type=str)
parser.add_argument("--rep", type=int)
args = parser.parse_args()

print("Loading variants")
variants = pd.read_csv(args.variants)
n_variants = len(variants)
print(f"Read {n_variants} variants")

print("Generating predictions")
os.environ["CUDA_VISIBLE_DEVICES"] = args.device
from decima.vep import predict_variant_effect

include_cols = [
    x for x in ["variant_id", "variant", "afc", "pip", "pos_variant"] if x in variants.columns
]

print(f"Replicate {args.rep}")
pq = os.path.join(args.out_dir, f"{args.prefix}_decima_rep{args.rep}.pq")
preds = predict_variant_effect(
        df_variant=variants,
        output_pq=pq,
        model=args.rep,
        device=0,
        gene_col="gene",
        include_cols=include_cols,
    )