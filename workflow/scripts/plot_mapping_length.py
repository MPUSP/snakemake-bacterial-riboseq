#!/usr/bin/python

# PLOT MAPPING LENGTH
# -----------------------------------------------------------------------------
#
# This script plots the length of all alignments contained in the bam file.
# Input bam file is checked for alignments.

import matplotlib
import pysam
import pandas as pd
from collections import Counter
from matplotlib import pyplot as plt


def get_aln_length(bamfile, sample):
    lengths = []
    for aln in bamfile.fetch():
        lengths.append(len(aln.query_sequence))
    bamfile.close()

    cnt_dict = Counter(lengths)
    df = pd.DataFrame.from_dict(cnt_dict, orient='index').sort_index()
    df.columns = [sample]
    return df


input = str(snakemake.input["bam"])
sample_name = snakemake.wildcards.sample
output_tsv = snakemake.output["tsv"]
output_pdf = snakemake.output["pdf"]
output_log = snakemake.log["path"]
log = []
error = []

bamfile = pysam.AlignmentFile(input, "rb")
num_aln = bamfile.count()

# check if bam file is empty
if num_aln != 0: 
    df = get_aln_length(bamfile=bamfile, sample=sample_name)

    # plot pdf
    matplotlib.rcParams['axes.grid'] = True

    x = df.index
    y = list(df[sample_name])
    zipped = zip(x, y)
    max_index, max_val = max(zipped, key=lambda x: x[1])

    fig, ax = plt.subplots()
    ax.plot(x, y, marker='o')
    ax.set_title('{}, mode at {}'.format(sample_name, max_index))
    ax.axvline(max_index, c='r', linestyle='--', alpha=0.3)
    ax.axvspan(20, 40, color='k', alpha=0.1)
    ax.set_ylabel('Counts')
    ax.set_xlabel('Read length')
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    plt.tight_layout()
    plt.savefig(output_pdf, bbox_inches='tight')

    # export count table
    df.index.set_names(['Length'], inplace=True)
    df.reset_index(inplace=True)
    df.to_csv(output_tsv, sep='\t', index=False)
    
    log += [f"{num_aln} alignments processed."]
else:
    # no alignment found.
    error += ["Bam file does not contain any alignment."]

# print error / log messages
if error:
    print("\n".join(error))
    raise ValueError(
        "Error occurred while processing BAM file, exit forced."
    )
else:
    log += ["Module finished successfully!\n"]
    log = ["PLOT_MAPPING_LENGTH: " + i for i in log]
    with open(output_log, "w") as log_file:
        log_file.write("\n".join(log))