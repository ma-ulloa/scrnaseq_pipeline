library(Seurat)
library(SoupX)
library(patchwork)
###############################################################
################## FUNCTIONS ##################################
###############################################################

# create function to read matrix, create seurat object


pp = function(filtered_path, sample_id){
  sobj = Read10X_h5(filtered_path)
  sobj = CreateSeuratObject(counts = sobj)
  sobj$sample_id         = sample_id
  sobj[["percent.mt"]]   = PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj[["percent.ribo"]] = PercentageFeatureSet(sobj, pattern = "^RP[SL]")
  return(sobj)
}


# Create basic clustering:
  
  
get_soup_groups = function(sobj){
  sobj = NormalizeData(sobj, verbose = FALSE)
  sobj = FindVariableFeatures(object = sobj,
                              nfeatures = 2000,
                              verbose = FALSE,
                              selection.method = "vst")
  sobj = ScaleData(sobj, verbose = FALSE)
  sobj = RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj = FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj = FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  sobj = RunUMAP(sobj, dims = 1:20, verbose = FALSE)
  return(sobj)
}


plot_soupgenes_umap = function(sobj_processed, soupgenes, plot_path){
  genes_to_plot = intersect(soupgenes, rownames(sobj_processed))
  if (length(genes_to_plot) == 0){
    cat("[WARN] No soupgenes found in the data — skipping UMAP plot\n")
    return(invisible(NULL))
  }

  plots    = FeaturePlot(sobj_processed, features = genes_to_plot, combine = FALSE)
  ncols    = min(4L, length(plots))
  nrows    = ceiling(length(plots) / ncols)
  combined = wrap_plots(plots, ncol = ncols)

  ggsave(plot_path, combined, width = 4 * ncols, height = 4 * nrows)
  cat("[INFO] Soupgene UMAP plot saved to:", plot_path, "\n")
}


add_soup_groups = function(sobj, soupgenes, plot_path){
  sobj_processed  = get_soup_groups(sobj)
  sobj$soup_group = sobj_processed@meta.data[["seurat_clusters"]]
  plot_soupgenes_umap(sobj_processed, soupgenes, plot_path)
  return(sobj)
}


make_soup = function(sobj, raw_path, soupgenes, fraction, apply_fractions, no_clusters){
  sample_id = as.character(sobj$sample_id[1])
  cat("=== ", sample_id, " ===\n\n")

  if (!(fraction %in% apply_fractions)){
    cat(sample_id, "— fraction '", fraction, "' not in apply_to_fractions, skipping SoupX correction\n")
    return(sobj)
  }

  if (!file.exists(raw_path))
    stop("Raw matrix not found: ", raw_path)

  raw = Read10X_h5(raw_path)
  sc  = SoupChannel(raw, sobj@assays$RNA$counts)
  sc  = setClusters(sc, sobj$soup_group)

  if (no_clusters){
    cat("Running without clusters (flagged in metadata)\n")
    useToEst = estimateNonExpressingCells(sc,
                                          nonExpressedGeneList = list(IG = soupgenes),
                                          clusters = FALSE)
  } else {
    useToEst = estimateNonExpressingCells(sc,
                                          nonExpressedGeneList = list(IG = soupgenes))
  }

  sc  = calculateContaminationFraction(sc, list(IG = soupgenes),
                                       useToEst = useToEst, forceAccept = TRUE)
  out = adjustCounts(sc)

  sobj[["original.counts"]] = CreateAssayObject(counts = sobj@assays$RNA$counts)
  sobj@assays$RNA$counts    = out
  cat("Done with ", sample_id, "\n")
  return(sobj)
}
