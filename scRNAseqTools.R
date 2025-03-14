# Useful functions for snRNAseq analysis

library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

# Get and plot counts ----
# Get counts
takeSoReturnCounts_spg <- function(so){

  cellTypeCounts <- so@meta.data %>% dplyr::count(sample, gt, group, subtype)

  }

# Plot counts 
plotCellTypeCounts_spg <- function(cellTypeCounts){

  cellTypeCounts %>% 
    ggplot(aes(x = group, y = n, color = gt)) +
    geom_boxplot(
      aes(color = gt), width = 0.5, size = 0.4,
      position = position_dodge(0.8)
    ) +
    geom_dotplot(
      aes(fill = gt, color = gt),
      binaxis='y', stackdir='center', dotsize = 0.8,
      position = position_dodge(0.8)
    ) + 
    facet_wrap(~scTypeClassificationSimple, scale = "free_y") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme_classic() +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.key.size = unit(1, "cm"))
  
}

# Generate deseq2 results ----
# https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

# Note this function is buggy and doesn't really work
pseudobulkReturnCountsAndMetadata_spg <- function(data, sampleGroupName = "sample", pseudobulkBySubGroup = TRUE, subGroupName = "subtype"){
  
  # browser()
  
  # Extract raw counts and metadata to create SingleCellExperiment object
  counts <- data[['RNA']]$counts
  
  metadata <- data@meta.data
  
  # Create single cell experiment object
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  
  if(pseudobulkBySubGroup){
    # Get cluster names (cell type names)
    cluster_names <- levels(as.factor(colData(sce)$subtype))
  }
  
  # Get sample names
  sample_names <- levels(as.factor(colData(sce)$sample))
  
  # Aggregate counts over samples
  if(pseudobulkBySubGroup){
   groups <- colData(sce)[, c(sampleGroupName, subGroupName)]
  }else{
    groups <- colData(sce)[, sampleGroupName]
  }
  
  # Aggregate across cluster-sample groups
  # transposing row/columns to have cell_ids as row names matching those of groups
  aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                  groupings = groups, fun = "sum") 
  
  # Split by cell types
  # Transpose aggregated matrix to have genes as rows and samples as columns
  aggr_counts <- t(aggr_counts)
  
  if(pseudobulkBySubGroup){
    # Loop over all cell types to extract corresponding counts, and store information in a list
    ## Initiate empty list
    counts_ls <- list()
    
    for (i in 1:length(cluster_names)) {
      
      ## Extract indexes of columns in the global matrix that match a given cluster
      column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[4]] == cluster_names[i])
      
      ## Store corresponding sub-matrix as one element of a list
      counts_ls[[i]] <- aggr_counts[, column_idx]
      names(counts_ls)[i] <- cluster_names[i]
      
    }
  }
  
  # Extract sample-level metadata
  metadata <- colData(sce) %>% 
    as.data.frame() %>% 
    dplyr::select(group, sample, gt)
  
  # Exclude duplicated rows
  metadata <- metadata[!duplicated(metadata), ]
  
  # Rename rows
  rownames(metadata) <- metadata$sample
  
  # Get number of cells per sample / cell type
  t <- table(colData(sce)$sample,
             colData(sce)$subtype)
  
  # Creating metadata list
  
  ## Initiate empty list
  metadata_ls <- list()
  
  for (i in 1:length(counts_ls)) {
    
    ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
    
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$subtype <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$sample  <- gsub("^[^_]*_", "", df$cluster_sample_id)
    
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df$cellType))
    cell_counts <- t[, idx]
    
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df$sample
                          , names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    
    ## Append cell_counts to data frame
    df$cell_count <- cell_counts
    
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    df <- plyr::join(df, metadata, 
                     by = intersect(names(df), names(metadata)))
    
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df) <- df$cluster_sample_id
    
    ## Store complete metadata for cluster i in list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df$cluster_id)
  
  }
  
  returnList <- list(counts = counts_ls, metadata = metadata_ls)
  
  return(returnList)
}



# Plot Deseq results ----
ggmaplotOnDeseq2Results_spg <- function(res, geneNames = NULL, top = 50, fdr = 0.05, fc = 2, size = 0.4){
  ggmaplot(res,
           fdr = fdr, fc = fc, size = size,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           genenames = geneNames,  
           legend = "top", top = top,
           font.label = c("bold", 11),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal()) + coord_flip()
}

pheatmapOfDeseq2Results_spg <- function(res, dds, 
                                        group = "WT_CFA|WT_D14", # Group is a string fed to grepL to pick the samples. These are in caps and can be combined with the or operator | 
                                        top = 100, log2fc = 1, padjusted = 0.05, genesToPlot = NULL, fontsize_row = 3, fontsize_col = 7){
  
  res.df <- as.data.frame(res)
  
  res.df$symbol <- rownames(res.df)
  
  if(!is.null(genesToPlot)){
    res.df.filtered <- res.df %>% filter(symbol %in% genesToPlot)}else{
    
    res.df.filtered <- res.df %>% filter(baseMean > 100, abs(log2FoldChange) > log2fc, padj < padjusted) %>% top_n(top, abs(log2FoldChange))
  }
  
  print(res.df.filtered)
  
  mat <- counts(dds, normalized = T)[rownames(res.df.filtered),] %>% 
    as.data.frame() %>% 
    select_if(grepl(group, colnames(.)))
  
  mat.z <- t(apply(mat, 1, scale))
  
  colnames(mat.z) <- colnames(mat)
  
  mat.z <- na.omit(mat.z)
  
  pheatmap(mat.z, fontsize_col = fontsize_col, fontsize_row = fontsize_row)
}

makeSortedGeneList_spg <- function(ddsResult) {
  df <- data.frame(ddsResult)
  df$gene <- rownames(ddsResult)
  df <- df %>% dplyr::select(gene, stat) %>% na.omit()
  df <- df[order(df$stat, decreasing = TRUE),]
  gl <- df[,2]
  names(gl) <- as.character(df[,1])
  gl <- sort(gl, decreasing = TRUE)
  return(gl)
}


# Sachin's stacked vln plot (ignores genes that are not in counts)
stackedVlnPlot_spg <- function (seurat, features, clus_ident = "seurat_clusters",
          assay = "RNA", layer = "counts", main = "", pt_size = 0) 
{
  Idents(seurat) <- clus_ident
  DefaultAssay(seurat) <- assay
  
  origNumberOfFeatures <- length(features)
  features <- features[features %in% rownames(seurat[["RNA"]])] # only keep features that are in the dataset
  
  features <- unique(features) # remove duplicates
  
  removedFeatureCount <- origNumberOfFeatures - length(features)
  if(removedFeatureCount > 0){print(paste0("Ignoring ", removedFeatureCount, " features that are not in the RNA assay"))}
  
  modify_vlnplot <- function(obj, feature, pt.size = pt_size, plot.margin = unit(c(-0.75, 
                                                                             0, -0.75, 0), "cm"), ...) {
    p <- VlnPlot(obj, features = feature, pt.size = pt.size, 
                 layer = layer, ...) + xlab("") + ylab(feature) + ggtitle("") + 
      theme(legend.position = "none", plot.title = element_blank(), 
            axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), axis.title.y = element_text(size = rel(1), 
                                                                        angle = 0), axis.text.y = element_text(size = rel(1)), 
            plot.margin = plot.margin)
    return(p)
  }
  extract_max <- function(p) {
    ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }
  Plot <- function(obj, features, pt.size = pt_size, plot.margin = unit(c(-0.75, 
                                                                    0, -0.75, 0), "cm"), ...) {
    plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, 
                                                                 feature = x, ...))
    plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                       vjust = 1), axis.ticks.x = element_line())
    ymaxs <- purrr::map_dbl(plot_list, extract_max)
    plot_list <- purrr::map2(plot_list, ymaxs, function(x, 
                                                        y) x + scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
    p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) + 
      patchwork::plot_annotation(title = main, theme = theme(plot.title = element_text(size = 18)))
    return(p)
  }
  return(Plot(obj = seurat, features = features))
}

# Make a function that runs sctransform, dim reduction on a seurat object
runSctransform_spg <- function(seurat, resolution = 0.05, dimNeighbor = 1:20, dimUmap = 1:20, verbose = TRUE){
  
  if(verbose){print("Running SCTransform")}
  
  seurat <- SCTransform(seurat)
  
  if(verbose){print("Running PCA")}
  
  seurat <- RunPCA(seurat)
  
  if(verbose){print("Running FindNeighbors and FindClusters")}
  
  seurat <- FindNeighbors(seurat, dims = dimNeighbor)
  
  seurat <- FindClusters(seurat, resolution = resolution)
  
  if(verbose){print("Running UMAP")}
  
  seurat <- RunUMAP(seurat, dims = dimUmap)
  
  return(seurat)
}

runPcaAndUmap_spg <- function(seurat, resolution = 0.05, dimNeighbor = 1:20, dimUmap = 1:20, verbose = TRUE){
  
  if(verbose){print("Running PCA")}
  
  seurat <- RunPCA(seurat)
  
  if(verbose){print("Running FindNeighbors and FindClusters")}
  
  seurat <- FindNeighbors(seurat, dims = dimNeighbor)
  
  clustree(seurat)
  
  seurat <- FindClusters(seurat, resolution = resolution)
  
  if(verbose){print("Running UMAP")}
  
  seurat <- RunUMAP(seurat, dims = dimUmap)
  
  return(seurat)
}

# function that subsets a pseudobulked SCE to a specific cluster then returns a DDS object. Will filter out any pseudobulked samples with less than 10 cells and warn.

ddsFromSCE_mds <- function(sce, clusterNames, clusterSelect, ddsDesign) {
  # Subset to just the clusters we want
  sce.select <- sce[,colData(sce)[[clusterNames]] %in% clusterSelect]
  colnames(sce.select) <- sce.select$sample
  # filter out any samples with less than 10 cells
  keep <- colData(sce.select)$ncells > 10
  # print a warning if dropping any samples
  if(sum(!keep) > 0) {
    message(paste0("Dropping ", sum(!keep), " samples with less than 10 cells for cluster ", clusterSelect))
  }
  # filter out the low cell samples
  sce.select <- sce.select[,keep]
  # Create a DESeqDataSet object
  dds <- DESeqDataSet(sce.select, design = ddsDesign)
  return(dds)
}

# function that makes a pca (pc1, pc2) of top 500 genes (default option) with DESeq2 based on given intgroup.
# it also makes a correleation heatmap ala DESeq2 vignette with pheatmap, using given variables for annotation.
# this is basically copy/paste of Sachin's exploratory code put into a function to make it easier to do repeatedly.
# returns a list of two plots, one is ggplot for pca, the other is pheatmap for correlation heatmap.
exploratoryDESeq2Analysis_mds <- function(dds, celltype, intGroup, hmapAnnotation) {
  dds
  rld <- rlog(dds, blind=TRUE)
  #pca <- DESeq2::plotPCA(rld, intgroup = intGroup) + ggtitle(celltype)
  pca <- pcaExplorer::pcaplot(rld, intgroup = intGroup) + ggtitle(celltype)
  hmap_cor <- pheatmap(cor(assay(rld)), annotation = as.data.frame(colData(rld))[,hmapAnnotation], drop=F, main = celltype, silent = TRUE)
  return(list(pca = pca, hmap_cor = hmap_cor))
}

exploratoryDESeq2Analysis_mds <- function(dds, celltype, intGroup, hmapAnnotation) {
  dds
  rld <- rlog(dds, blind=TRUE)
  pca <- pcaExplorer::pcaplot(rld, intgroup = intGroup) + ggtitle(celltype)
  hmap_cor <- pheatmap(cor(assay(rld)), annotation = as.data.frame(colData(rld)[, hmapAnnotation, drop=FALSE]), main = celltype, silent = TRUE)
  return(list(pca = pca, hmap_cor = hmap_cor))
}

# function that extracts from a Seurat data the metadata and raw rna counts then stores them in a new seurat object. basically what DietSeurat should do but doesn't.

cleanSeurat <- function(so) {
  md <- so@meta.data
  new <- CreateSeuratObject(counts = LayerData(so, assay = "RNA", layer = "counts"))
  new@meta.data <- md
  return(new)
}

# function to make a stacked barplot off of a Seruat metadata object. returns ggplot object
## to make this more broadly applicable, I should add some arguments on how to do colors. right now its hard coded scCustomize polychrome which is 36 or less
stackedBarPlot_mds <- function(df, x, y, order_x=NULL, order_y=NULL, title) {
  # Calculate proportions
  df_proportions <- df %>%
    group_by({{ x }}, {{ y }}) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by({{ x }}) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # Order the levels of the x variable
  if(!is.null(order_x)) {
    df_proportions <- df_proportions %>%
      mutate({{ x }} := factor({{ x }}, levels = order_x))
  }
  
  # order the levels of the y variable
  if(!(is.null(order_y))) {
    df_proportions <- df_proportions %>%
      mutate({{ y }} := factor({{ y }}, levels = order_y))
  }
  
  
  # Determine the number of levels in y
  num_colors <- length(unique(df_proportions[[deparse(substitute(y))]]))  # Count unique levels
  
  # Make the ggplot
  p <- ggplot(df_proportions, aes(x = {{ x }}, y = proportion, fill = {{ y }})) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent_format()) +  # Format y-axis as percentages
    labs(x = rlang::as_string(quo_name(enquo(x))), y = "Proportion", fill = rlang::as_string(quo_name(enquo(y)))) +  # Convert quosure to string
    theme_minimal() +
    scale_fill_manual(values = DiscretePalette_scCustomize(num_colors = num_colors, palette = "polychrome")) +  # Use the custom colors
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.text = element_text(size = 14)) + # cosmetics - rotate labels, text size
    ggtitle(title)
  
  return(p)
}

# function that extracts named gene lists out of ego object
extractGeneListsFromEgo_mds <- function(ego) {
  go_genes <- ego@result$geneID # grab genes from geneID column
  go_genes_list <- strsplit(go_genes, "/") # split the genes by go term by / character
  names(go_genes_list) <- ego@result$Description # name the gene sets
  return(go_genes_list)
}

# function that takes a df from deseq2 results and an ego object and makes a cnetplot with genes colored by fold change
# the default name for foldChangeName is the one DESeq2 uses, edgeR uses LogFC
cnetplotFromDf_mds <- function(df, ego, pathways, foldChangeName = "log2FoldChange") {
  log2fc <- df[[foldChangeName]]
  names(log2fc) <- rownames(df)
  p <- enrichplot::cnetplot(ego, color.params = list(foldChange = log2fc), showCategory = pathways)
  return(p)
}

# function that takes a dataframe generated by topTags (edgeR) and a list of gene sets and makes a limma::barcode plot using signed F score as x-axis
# it optionally can show log2foldchange weights (1 gene set) or 2 gene sets with the arguments. no return as limma::barcodeplot plot is a side effect, not a value
makeBarcodePlot_mds <- function(res, ids, title1, title2 = NULL, weights = NULL, ...) {
  res <- na.omit(res)
  res$signedF <- ifelse(res$logFC > 0, res$F, -res$F)
  idx <- ids2indices(ids, rownames(res))
  weights.label = NULL
  if(!is.null(weights)) {
    weights = res$logFC[idx[[title1]]]
    weights.label = "log2FC"
  }
  if(is.null(title2)) {
    limma::barcodeplot(statistics = res$signedF, index = idx[[title1]], gene.weights = weights, weights.label = weights.label, xlab = "Signed F Statistic", ...)
  }
  else {
    limma::barcodeplot(statistics = res$signedF, index = idx[[title1]], index2 = idx[[title2]], gene.weights = weights, xlab = "Signed F Statistic", ...)
  }
}

# function that runs camera on a results object returned by scran::pseudoBulkDGE. It's really running it on the returned DGEGLM.
# CAMERA automatically uses last column of design matrix unless you specify another matrix, which has to be a numeric, can't use characters to feed into makeContrasts
runCamera_mds <- function(res, ids, contrast) {
  fit <- metadata(res)$fit
  idx <- ids2indices(ids, rownames(fit))
  x <- camera(y = fit, index = idx, contrast = contrast)
  return(x)
}

 
# function that takes a list of df: res, up, down. it then runs clusterProfiler::enrichGO on the up and down gene lists and returns a list of results. It checks to make sure there are actual genes in each df, and always returns a list
# i specify the universe in the function call isntead of by grabbing the rownames of res because we filter out NAs from the res object? Maybe we shouldn't do that? But leaving it as is now.
runEnrichGoOnList_mds <- function(df.ls, ontology = "BP", universe = universe, OrgDb = org.Mm.eg.db, ...) {

  # Initialize ego.up and ego.down as NULL. This is so we can check if they're still null at the end.
  ego.up <- NULL
  ego.down <- NULL
  
  # check ot see if ther eare any genes in the up gene list
  if (nrow(df.ls$up) > 0) {
    # run enrichGO on the up gene list
    ego.up <- enrichGO(rownames(df.ls$up), keyType = "SYMBOL", OrgDb = OrgDb, ont = ontology, universe = universe, ...)
    # if there were any enriched pathways, run pairwise_termsim on the ego object
    if (dim(ego.up)[1] > 0) {
      ego.up <- enrichplot::pairwise_termsim(ego.up)
    }
  }
  # check to see if ther eare any genes in the down gene list
  if (nrow(df.ls$down) > 0) {
    # run enrichGO on the down gene list
    ego.down <- enrichGO(rownames(df.ls$down), keyType = "SYMBOL", OrgDb = OrgDb, ont = ontology, universe = universe, ...)
    # if there are any enriched pathways, run pairwise_termsim
    if (dim(ego.down)[1] > 0) {
      ego.down <- enrichplot::pairwise_termsim(ego.down)
    }
  }
  
  # Return a list with ego.up and ego.down only if they were created
  result <- list()
  if (!is.null(ego.up)) result$ego.up <- ego.up
  if (!is.null(ego.down)) result$ego.down <- ego.down
  return(result)
}

# write a function that runs clusterProfiler::enricher, this works exactly the same as as the enrichGo function above, but now we need to provide a term2gene. we also have to make sure at least 1 gene maps to term2gene otherwise it will throw an error.
runEnricherOnLists_mds <- function(res, term2gene, universe, ...) {
  enrich.up <- NULL
  enrich.down <- NULL
  
  if(nrow(res$up) > 0) {
    # check ot make sure some of the gnes are int he term2gene, otherwise it will throw an error about being unable to map
    if(sum(rownames(res$up) %in% term2gene$gene) > 0) {
      enrich.up <- enricher(gene = rownames(res$up), universe = universe, TERM2GENE = term2gene, ...)
      if(dim(enrich.up)[1] > 0) {
        enrich.up <- pairwise_termsim(enrich.up)
      }
    }
  }
  if(nrow(res$down) > 0) {
    # check ot make sure some of the gnes are int he term2gene, otherwise it will throw an error about being unable to map
    if(sum(rownames(res$down) %in% term2gene$gene) > 0) {
      enrich.down <- enricher(gene = rownames(res$down), universe = universe, TERM2GENE = term2gene, ...)
      if(dim(enrich.down)[1] > 0) {
        enrich.down <- pairwise_termsim(enrich.down)
      }
    }
  }
  result <- list()
  if(!is.null(enrich.up)) result$enrich.up <- enrich.up
  if(!is.null(enrich.down)) result$enrich.down <- enrich.down
  return(result)
}

# function that takes a list list of bulk RNAseq results where rows are gene names, and then returns a list with all of them (after na.omit), then just the up and just the down
subsetDegList_mds <- function(result, foldChangeName = "logFC", foldChangeCutoff = (1.5), pvalName = "FDR", pvalCutoff = 0.05, expressionCutoffName = "logCPM", expressionCutoff = NULL) {
  # subset the results to only include genes that pass the fold change and p value cutoffs
  x <- na.omit(result) #both DESeq2 and edgeR downstream of pseudoBulkDGE return resutls with NAs in them, so get rid of those.
  cutoff <- log2(foldChangeCutoff)
  # if an expression cutoff is set, use it, otherwise don't filter on expression and just on logFC and pvalue
  if(!is.null(expressionCutoff)) {
    up <- x[x[[foldChangeName]] > cutoff & x[[pvalName]] < pvalCutoff & x[[expressionCutoffName]] > expressionCutoff,]
    down <- x[x[[foldChangeName]] < -cutoff & x[[pvalName]] < pvalCutoff & x[[expressionCutoffName]] > expressionCutoff,]
  } else {
    up <- x[x[[foldChangeName]] > cutoff & x[[pvalName]] < pvalCutoff,]
    down <- x[x[[foldChangeName]] < -cutoff & x[[pvalName]] < pvalCutoff,]
  }
  return(list(res = x, up = up, down = down))
}

# function that takes a gene list, a df, and something to wory by. it subsetes the df to match the gene list, sorts it by the specified argument, and returns that sorted subsetted list of genes.
subset_and_sort_genes <- function(gene_list, df, sort_by = "logFC") {
  # subset to just the genes
  subset_df <- df[rownames(df) %in% gene_list, ]
  # Sort the subset dataframe by something in descending order
  sorted_df <- subset_df[order(-subset_df[[sort_by]]), ]
  # Return the sorted list of gene names
  return(rownames(sorted_df))
}

# function that takes a DataFrame generated by edgeR (pseudobulkDGE) and returns the PValue where FDR is greater than threshold. Could work on test statistic as well but woudl have to make sure signs ok
getPValueFromFDR_mds <- function(DF, pValueName = "PValue", FDRName = "FDR", FDRcutoff = 0.05, decreasing = FALSE) {
  df <- as.data.frame(DF)
  #df <- df[order(df[[FDRName]], decreasing = decreasing),]
  first_non_significant_PValue <- df %>%
    dplyr::filter(.data[[FDRName]] > 0.05) %>% #.data to force FDRName to be variable holding column name
    dplyr::slice(1) %>%
    dplyr::pull(pValueName)
  return(first_non_significant_PValue)
}

calcPercCellsExpressing_mds <- function(gene, so) {
  total_cells <- ncol(so)
  cells_expressing_gene <- sum(GetAssayData(soNeuron.no42, layer = "counts", assay = "RNA")[gene, ] > 0)
  percent_cells_expressing_gene <- (cells_expressing_gene / total_cells) * 100
  return(percent_cells_expressing_gene)
}

subsetGenesToPercentExpressed_mds <- function(so, genes, threshold) {
  # get counts
  counts <- GetAssayData(so, assay = "RNA", layer = "counts")
  # filter counts to input genes
  counts <- counts[genes,]
  # calculate percent expression
  genes.percent.expression <- rowMeans(counts > 0) * 100
  # filter counts to just the genes that pass the threshold
  counts <- counts[genes.percent.expression > threshold,]
  # order by percent expressed
  filtered_genes <- rownames(counts)
  ordered_genes <- filtered_genes[order(genes.percent.expression[filtered_genes], decreasing = TRUE)]
  # return the filtered genes
  return(ordered_genes)
}


# function for calculating distances and making heatmap, optionally allows assignign weights to control how the dendrogram arms swing. data should have samples in columns and genes on row (like all tx)
# assumes you want to cluster samples because otherwise why use this function. also allows for transposing so the heatmap has the genes on columns and samples on rows, but the data has to be put in like normal
# extension of the function I wrote for MS187 nlrx1ko astrocyte in vitro bulk rnaseq
# extra arguments go to pheatmap
doHeatmap_mds <- function(data, genes, weights = NULL, distanceAlg = "euclidean", linkageMethod = "complete", clustGenes = TRUE, transpose = FALSE, ...) {
  # subset input data to just the genes of interest
  data <- data[genes,]
  # determine if distance algorithm is spearman or pearson, in which case it's 1-. And also, you transpose data for rows, not columns.
  if(distanceAlg %in% c("spearman", "pearson")) {
    dist_genes <- as.dist(1-cor(t(data), method = distanceAlg))
    dist_samples <- as.dist(1-cor(data, method = distanceAlg))
    # if you're using euclidean or manhattan, you don't do 1-, and you transpose for columns
  } else {
    dist_genes <- dist(data, method = distanceAlg)
    dist_samples <- dist(t(data), method = distanceAlg)
  }
  # with the distance calculations, then do heirarchical clustering using specified linkage method
  clust_genes <- hclust(dist_genes, method = linkageMethod)
  clust_samples <- hclust(dist_samples, method = linkageMethod)
  # to use weights, for columns convert hclust to dendrogram, then reorder it.
  dend_samples <- as.dendrogram(clust_samples)
  if(!is.null(weights)) {
    dend_samples <- reorder(dend_samples, wts = weights, agglo.FUN = "mean")
  }

  # plot the heatmap with pheatmap. Scale for gene optionally, have to convert the column dendrogram back to hclust.
  # this changes based on some of the options
  if(clustGenes == TRUE) {
    if(transpose == FALSE) {
      pheatmap(as.matrix(data),
               cluster_rows = as.hclust(clust_genes),
               cluster_cols = as.hclust(dend_samples), scale="row", ...)
    } else {
      pheatmap(as.matrix(t(data)),
               cluster_cols = as.hclust(clust_genes),
               cluster_rows = as.hclust(dend_samples), scale="column", ...)
    }
  }
  else {
    if(transpose == FALSE) {
      pheatmap(as.matrix(data),
               cluster_rows = FALSE,
               cluster_cols = as.hclust(dend_samples), scale="row", ...)
    } else {
      pheatmap(as.matrix(t(data)),
               cluster_cols = FALSE,
               cluster_rows = as.hclust(dend_samples), scale="column", ...)
    
    }
  }
}
