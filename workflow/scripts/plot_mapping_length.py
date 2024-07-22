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


def get_aln_length(bamfile, outfolder, sample):
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

    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-b', dest='bam', required=True, help='bam file')
    parser.add_argument('-o', dest='output_folder', required=True, help='output folder')
    parser.add_argument('-s', dest='sample_name', required=True, help='sample name')

    args = vars(parser.parse_args())
    file = args['bam']
    sample_name = args['sample_name']
    outfolder=args['output_folder']


    try:
        if not file.endswith('bam'):
            sys.stderr.write("Input file needs to be in bam format. " +
                             "Exit forced\n.")
            sys.exit(1)
        elif not os.path.exists(file):
            sys.stderr.write("Input file %s does not exist. " +
                             "Exit forced\n." % file)
            sys.exit(1)
        bamfile = pysam.AlignmentFile(file, "rb")
        df = get_aln_length(bamfile=bamfile, outfolder=outfolder,
                            sample=sample_name)

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
        plt.savefig(os.path.join(outfolder, f'{sample_name}_mapped_length.pdf'), bbox_inches='tight')

        # export coount table
        df.index.set_names(['Length'], inplace=True)
        df.reset_index(inplace=True)
        df.to_csv(os.path.join(outfolder, f'{sample_name}_length_dist.tsv'), sep='\t', index=False)

    except:
        sys.stderr.write("Error ocurred when processing bam file: %s.\n"
                         % file)
        raise RuntimeError
        sys.exit(1)
