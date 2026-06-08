#!/usr/bin/env python3
"""
02_multisample_report.py
────────────────────────────────────────────────────────────────
Joins qc metrics of multiple samples into one report

Loads data, merges QC metrics tables, builds figures, and writes
an HTML + PDF report. Designed to run pre- and post-filtering.

Usage:
    python 02_multisample_report.py \
        --input      \
        --metadata  config/metadata.csv \
        --config    config/config.yaml \
        --sample_id KS2103T1_1 \
        --stage     pre \
        --output    results/qc/pre_KS2103T1_1_qc.html

Input format is auto-detected from the file extension:
    .h5ad → AnnData  |  .h5 → 10x or generic HDF5
    .mtx / dir → 10x MEX  |  .csv / .txt / .tsv → dense matrix
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import scanpy as sc
import yaml

# Pipeline modules