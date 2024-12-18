import numpy as np
import os
import modiscolite
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-seq_file", 
                    help="Path to the one hot sequences")
parser.add_argument("-attr_file", 
                    help="Path to attributions")
parser.add_argument("-meme_file", 
                    help="Where the meme motifs for scanning are stored")
parser.add_argument("-out_dir", 
                    help="Directory where to save results")
args = parser.parse_args()

seq_file = args.seq_file
attr_file = args.attr_file
meme_file = args.meme_file
out_dir = args.out_dir

print("Loading data")
sequences = np.load(seq_file)
attributions = np.load(attr_file)

print("Running modisco")
sequences = sequences.transpose(0, 2, 1).astype("float32")
attributions = attributions.transpose(0, 2, 1).astype("float32")
pos_patterns, neg_patterns = modiscolite.tfmodisco.TFMoDISco(
    hypothetical_contribs=attributions,
    one_hot=sequences,
    max_seqlets_per_metacluster=10000,
)

print("Writing modisco output")
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