#!/usr/bin/python

# EXTRACT FEATURES 
# -----------------------------------------------------------------------------
#
# This script extracts selected biotype features from input GFF file.


from os import path
from BCBio.GFF import GFFExaminer
from BCBio import GFF

input_gff = snakemake.input["gff"]
output_gff = snakemake.output["gff"]
feature_list = snakemake.params["features"]
output_log = snakemake.log["path"]
log = []
error = []

if not path.exists(input_gff):
        error += ["The parameter 'gff' is not a valid path to a GFF file"]

try:
    with open(input_gff, 'r') as gff_file:
        examiner = GFFExaminer()
        limits = dict(gff_source_type=[('RefSeq', 'gene')])
        for rec in GFF.parse(gff_file, limit_info=limits):
            sel_features = []
            for feat in rec.features:
                if 'gene_biotype' in feat.qualifiers.keys():
                    if feat.qualifiers['gene_biotype'][0] in feature_list:    
                        sel_features += [feat]
            rec.features = sel_features
            if len(sel_features) == 0:
                error += ["No features found!"]
except IOError:
    error += [f"Supplied GFF file '{input_gff}' can not be opened"]

try:
    with open(output_gff, 'w') as gff_out:
        GFF.write([rec], gff_out)
except (IOError, ValueError, EOFError) as e:
    error += [e]

# print error/log messages
if error:
    print("\n".join(error))
    raise ValueError(
        "Either location of supplied GFF file was not correct or list of features were not found, quitting"
    )
else:
    log += ["Module finished successfully\n"]
    log = ["EXTRACT FEATURES: " + i for i in log]
    with open(output_log, "w") as log_file:
        log_file.write("\n".join(log))