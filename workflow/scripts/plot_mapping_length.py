#!/usr/bin/python

import os
import sys
import argparse
import matplotlib
import pysam
import pandas as pd
from collections import Counter
from matplotlib import pyplot as plt


__authors__ = "Rina Ahmed-Begrich"
__mail__ = "begrich@mpusp.mpg.de"
__program_name__ = "plot_mapping_length.py"
__version_info__ = ('0', '0', '1')
__version__ = '.'.join(__version_info__)


def get_aln_length(bamfile, sample):
    lengths = []
    for aln in bamfile.fetch():
        lengths.append(len(aln.query_sequence))
    bamfile.close()

    cnt_dict = Counter(lengths)
    df = pd.DataFrame.from_dict(cnt_dict, orient='index').sort_index()
    df.columns=[sample]
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="""
Example Usage:

    plot_mapping_length.py -b <bam_file> -s <sample_name> -o <output_folder>

Written by %s (%s),
Max Planck Unit for the Science of Pathogens. Copyright (c) 2024.
Copyright Holder All Rights Reserved.
""" % (__authors__, __mail__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    input = str(snakemake.input["bam"])
    sample_name = snakemake.wildcards.sample
    output_tsv = snakemake.output["tsv"]
    output_pdf = snakemake.output["pdf"]

    try:
        if not os.path.exists(input):
            sys.stderr.write(
                "Input file %s does not exist. " + "Exit forced\n." % input
            )
            sys.exit(1)
        bamfile = pysam.AlignmentFile(input, "rb")
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

        # export coount table
        df.index.set_names(['Length'], inplace=True)
        df.reset_index(inplace=True)
        df.to_csv(output_tsv, sep='\t', index=False)

    except:
        sys.stderr.write("Error ocurred when processing bam file: %s.\n"
                         % input)
        raise RuntimeError
        sys.exit(1)
