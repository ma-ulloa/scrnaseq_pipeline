library(here)
library(SoupX)
library(Seurat)
library(anndataR)
library(tidyverse)
library(argparse)
source(here("workflow/scripts/utils/soupx_utils.R"))


###############################################################
################## ARGUMENTS ##################################
###############################################################

parser = ArgumentParser(description = "Ambient RNA correction with SoupX — per sample")
parser$add_argument("--filtered",   required = TRUE,
                    help = "Path to filtered_feature_bc_matrix.h5")
parser$add_argument("--raw",        required = TRUE,
                    help = "Path to raw_feature_bc_matrix.h5")
parser$add_argument("--soupgenes",  required = TRUE,
                    help = "Path to text file with one soup marker gene per line")
parser$add_argument("--metadata",   required = TRUE,
                    help = "Path to metadata CSV (first column = sample_id)")
parser$add_argument("--sample_id",  required = TRUE,
                    help = "Sample ID to look up in metadata")
parser$add_argument("--output",     required = TRUE,
                    help = "Output .h5ad path")
parser$add_argument("--apply_to_cellss", required = TRUE, nargs = "+",
                    help = "cells(s) to apply SoupX to (e.g. Immune All). Samples not in this list are returned uncorrected.")
args = parser$parse_args()

soupgenes  = readLines(args$soupgenes)
metadata   = read.csv(args$metadata, stringsAsFactors = FALSE)
sample_row = metadata[metadata[[1]] == args$sample_id, ]
if (nrow(sample_row) == 0) stop("sample_id '", args$sample_id, "' not found in metadata")
cells    = sample_row$cells
no_clusters = as.logical(sample_row$no_clusters)


###############################################################
################## MAIN #######################################
###############################################################

cat("[INFO] Loading filtered matrix:", args$filtered, "\n")
sobj = pp(args$filtered, args$sample_id)

cat("[INFO] Computing clustering for SoupX\n")
umap_pre_path = sub("\\.h5ad$", "_soupgenes_umap_pre.pdf", args$output)
sobj = add_soup_groups(sobj, soupgenes, umap_pre_path)

cat("[INFO] Running SoupX correction\n")
sobj = make_soup(sobj, args$raw, soupgenes, cells, args$apply_to_cellss, no_clusters)

cat("[INFO] Plotting post-correction soupgene UMAP\n")
umap_post_path  = sub("\\.h5ad$", "_soupgenes_umap_post.pdf", args$output)
sobj_post       = get_soup_groups(sobj)
plot_soupgenes_umap(sobj_post, soupgenes, umap_post_path)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
sobj_sce = as.SingleCellExperiment(sobj)
anndataR::write_h5ad(sobj_sce, args$output)
cat("[INFO] Saved:", args$output, "\n")
