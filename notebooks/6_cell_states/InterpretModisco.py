import numpy as np
import pandas as pd
import os
import modiscolite

# Parameters

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-window", 
                    help="Window around tss to scan/modisco",
                    type=int, default=10000)
parser.add_argument("-seq_file", 
                    help="Path to the one hot sequences")
parser.add_argument("-attr_file", 
                    help="Path to attributions")
parser.add_argument("-meme_file", 
                    help="Where the meme motifs for scanning are stored")
parser.add_argument("-out_dir", 
                    help="Directory where to save results")
parser.add_argument("-seqlets", 
                    help="How many seqlets to use", type=int, default=10_000)
args = parser.parse_args()

seq_file = args.seq_file
attr_file = args.attr_file
meme_file = args.meme_file
out_dir = args.out_dir
window = args.window
seqlets = args.seqlets

seq_len = 524288
tss_pos = 163840

window_start = tss_pos - window
window_end = tss_pos + window

# Load and Prepare data
sequences = np.load(seq_file)
attributions = np.load(attr_file)

print("Modisco-ing")

# Madjdandzic et al. trick
attributions = attributions - attributions.mean(1,keepdims=True)

# Trim attributions
sequences = sequences[:,:,window_start:window_end]
attributions = attributions[:,:,window_start:window_end]

print("Running modisco")
one_hot_arr = sequences.transpose(0, 2, 1).astype("float32")
attrs = attributions.transpose(0, 2, 1).astype("float32")
pos_patterns, neg_patterns = modiscolite.tfmodisco.TFMoDISco(
    hypothetical_contribs=attrs, one_hot=one_hot_arr, max_seqlets_per_metacluster=seqlets)

print("Writing modisco output")
out_dir = os.path.join(out_dir,'modisco_full')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

h5_file = os.path.join(out_dir, "modisco_report.h5")
modiscolite.io.save_hdf5(h5_file, pos_patterns, neg_patterns, window_size=20)

print("Making report")
modiscolite.report.report_motifs(
    h5_file,
    out_dir,
    is_writing_tomtom_matrix=True,
    top_n_matches=10,
    meme_motif_db=meme_file,
    img_path_suffix="./",
    trim_threshold=0.2,
)

print("Done")