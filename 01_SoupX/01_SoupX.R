library(here)
library(SoupX)
library(Seurat)
library(tidyverse)
library(readxl)
data_path = here("data")
matrices_folders = list.dirs(here(data_path, "matrices"),
                              full.names = T)[-1]

###############################################################
################## FUNCTIONS ##################################
###############################################################

# create function to read matrix, create seurat object


pp = function(matrix_folder){
  path = file.path(matrix_folder,"filtered_feature_bc_matrix.h5")
  sobj = Read10X_h5(path)
  sobj = CreateSeuratObject(counts = sobj)
  sobj$sample_id = gsub(pattern = "_features",
                        replacement = "",
                        x = basename(matrix_folder))
  # add QC
  sobj[["percent.mt"]] = PercentageFeatureSet(sobj, pattern = "^MT-")
  
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
  sobj = ScaleData(sobj , verbose = F)
  sobj = RunPCA(sobj, npcs = 20, verbose = F)
  sobj = FindNeighbors(sobj, dims = 1:20, verbose = F)
  sobj = FindClusters(sobj, resolution = 0.5, verbose = F)
  
  return(sobj@meta.data[["seurat_clusters"]])
}


add_soup_groups = function(sobj){
  sobj$soup_group = get_soup_groups(sobj)
  return(sobj)
}




# Make soup:
  
# c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            # "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")

soupgenes = snakemake@input$soupgenes 




make_soup = function(sobj){
  
  sample_id = as.character(sobj$sample_id[1])
  cat("=== ", sample_id, " ===\n\n")
  if(phenoData[phenoData$sample_name == unique(sobj$sample_id),"cells"] == "EPI+DN"){
    cat(sample_id," is EPIDN\n")
  }
  
  else{
    
    path = file.path(data_path,"matrices",paste0(unique(sobj$sample_id),"_features"),"raw_feature_bc_matrix.h5")
    # Ensure path is a single string
    if(length(path) != 1 || !file.exists(path)){
      stop("Path does not exist or is not a single string: ", path,"\n")
    }
    
    raw = Read10X_h5(path)
    
    sc = SoupChannel(raw, sobj@assays$RNA$counts)
    sc = setClusters(sc, sobj$soup_group)
    
    if(unique(sobj$sample_id) == "KS2204T1_2"){
      cat("Running sample KS2204T1_2 without clusters because high contamination \n")
      useToEst = estimateNonExpressingCells(sc,
                                            nonExpressedGeneList = list(IG = igGenes),
                                            clusters = F)
    }

   else{
      useToEst = estimateNonExpressingCells(sc,
                                            nonExpressedGeneList = list(IG = igGenes))
    }
    
    
    sc = calculateContaminationFraction(sc, 
                                        list(IG = igGenes), 
                                        useToEst = useToEst,
                                        forceAccept = TRUE)
    out = adjustCounts(sc)
    
    sobj[["original.counts"]] = CreateAssayObject(counts = sobj@assays$RNA$counts)
    sobj@assays$RNA$counts <- out
    cat("Done with ", sample_id, "\n")
  }
  
  return(sobj)
  
}  










# create vector of seurat objects

data_list = sapply(matrices_folders, pp)

names(data_list) = gsub(pattern = "_features",
                        replacement = "",
                        x = basename(names(data_list)))

# remove descending 

descending_names = phenoData$sample_name[phenoData$region == "descendent_colon"]

data_list = data_list[names(data_list) != descending_names]
data_list = sapply(data_list , add_soup_groups)

data_list = sapply(data_list,make_soup)



# save
saveRDS(data_list, here("results/0_soupX/data_list.rds"))


