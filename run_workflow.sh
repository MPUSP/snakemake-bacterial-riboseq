#!/bin/bash

THREADS=40

snakemake --cores ${THREADS} --conda-frontend mamba --use-conda --rerun-incomplet --conda-cleanup-pkgs cache --printshellcmds
