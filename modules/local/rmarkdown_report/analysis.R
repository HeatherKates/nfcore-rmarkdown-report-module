# Load required libraries
library(limma)
library(edgeR)
library(tidyverse)
library(openxlsx)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(RColorBrewer)
library(ggrepel)
library(hues)
library(Biobase)
library(stringr)
library(ggfortify)
library(qs)
library(DT)
library(dplyr)
library(NOISeq)
library(AnnotationDbi)
library(org.Mm.eg.db)
source("HRK_funcs.R")
select <- dplyr::select
annotation_obj <- get(report_params$annotation_db, envir = asNamespace(report_params$annotation_db)) 

tmp_out_dir <- tempfile("rmd_tmpdir_")
dir.create(tmp_out_dir)

# Save it in the list of parameters
report_params$out_dir <- tmp_out_dir

# Step 1: List RSEM files
files.list <- list.files(path = report_params$rsem_dir, pattern = 'genes.results$', full.names = TRUE)

# Step 2: Extract sample names from filenames (remove suffix)
file_sample_names <- gsub(".genes.results$", "", basename(files.list))

# Step 3: Create data frame mapping sample names to files
file_df <- data.frame(SampleName = file_sample_names, File = files.list, stringsAsFactors = FALSE)

# Step 4: Load sample keys
sample.keys <- read_csv(report_params$sample_data)

# Step 5: Filter and re-order file_df to match sample.keys
file_df <- file_df[file_df$SampleName %in% sample.keys$SampleName, ]
file_df <- file_df[match(sample.keys$SampleName, file_df$SampleName), ]  # ensures correct order

# Step 6: Read only matched files
DGE <- readDGE(files = file_df$File, columns = c(1, 5))

# Step 7: Rename samples
colnames(DGE$counts) <- file_df$SampleName
rownames(DGE$samples) <- file_df$SampleName
DGE$samples$SampleName <- file_df$SampleName

# Step 8: Merge metadata (safe way that preserves order)
DGE$samples <- left_join(DGE$samples, sample.keys, by = "SampleName")
rownames(DGE$samples) <- DGE$samples$SampleName

# Annotate genes
ensemblid <- rownames(DGE)

# Get gene annotations
genes <- AnnotationDbi::select(annotation_obj, 
                               keys = ensemblid, 
                               columns = c("ENSEMBL", "SYMBOL"), 
                               keytype = "ENSEMBL")

# Remove duplicated Ensembl IDs (keeping the first occurrence)
genes <- genes[!duplicated(genes$ENSEMBL), ]

# Ensure we retain all IDs and row order in DGE
matched_genes <- genes[match(ensemblid, genes$ENSEMBL), ]

# Assign the matched gene annotations as a new column in DGE
DGE$genes <- matched_genes

# Define treatments
treatment.all <- as.factor(DGE$samples[[report_params$group_var]])

# Filter genes based on expression level (edgeR)
keep.exprs <- filterByExpr(DGE, group = treatment.all)
DGE.edgeRfilt <- DGE[keep.exprs, , keep.lib.sizes = FALSE] 

# Normalize using TMM method (edgeR)
DGE.edgeRfilt <- calcNormFactors(DGE.edgeRfilt, method = 'TMM')

# Filter genes based on NOISeq
DGE.NOIseqfilt <- DGE
DGE.NOIseqfilt$counts <- filtered.data(
  DGE$counts, 
  factor = DGE$samples[[report_params$group_var]], 
  norm = FALSE, 
  depth = NULL, 
  method = 1, 
  cv.cutoff = 100, 
  cpm = 1, 
  p.adj = "fdr"
)

DGE.NOIseqfilt$samples$lib.size <- apply(DGE.NOIseqfilt$counts, 2, sum)
DGE.NOIseqfilt <- calcNormFactors(DGE.NOIseqfilt, method = 'TMM')

# Create design matriDGE
design.mat <- model.matrix(~ 0 + DGE$samples[[report_params$group_var]])
rownames(design.mat) <- DGE$samples$SampleName
colnames(design.mat) <- make.names(levels(as.factor(DGE$samples[[report_params$group_var]])))

# Perform voom transformation
# Create voom plot and save voom object

png(filename = paste0(report_params$out_dir, "/voom_plot.png"))  # Open PNG device
v <- voom(DGE.NOIseqfilt, design.mat, plot = TRUE)
dev.off()  # Close PNG device

# Fit linear model
vfit <- lmFit(v, design.mat)

# Extract unique contrast names
# contrast_strings <- unique(unlist(strsplit(DGE.NOIseqfilt$samples$Contrast, ";")))
# 3. Extract contrasts list
if (!is.null(report_params$contrasts) && file.exists(report_params$contrasts)) {
  contrast_strings <- readLines(report_params$contrasts)
  contrast_strings <- contrast_strings[!grepl("^\\s*#", contrast_strings)]
  contrast_strings <- contrast_strings[contrast_strings != ""]
  contrast_strings <- str_replace_all(contrast_strings, "_vs_", "-")
} else {
  contrast_strings <- unique(unlist(strsplit(sample.keys$Contrast, ";")))
}


# Generate contrast list
contrasts_list <- setNames(
  lapply(contrast_strings, function(contrast) {
    groups <- trimws(unlist(strsplit(contrast, "-")))  # Remove whitespace
    
    print(paste("Processing contrast:", paste(groups, collapse = " - ")))  # Debugging step
    
    # Directly pass the contrast formula (no eval/parse needed!)
    makeContrasts(contrasts = paste(groups[1], "-", groups[2]), 
                  levels = colnames(vfit$design))
  }),
  paste0(gsub("-", "_vs_", contrast_strings))  # Clean names for list elements
)


# Print the list to check
print(contrasts_list)


# Compute contrasts and apply eBayes
efit_list <- lapply(contrasts_list, function(contr) {
  vfit_contr <- contrasts.fit(vfit, contr)
  eBayes(vfit_contr)
})

# Save differential expression results

library(clusterProfiler)

# Create result data frames dynamically
efit_results_list <- setNames(
  lapply(seq_along(efit_list), function(i) {
    create_efit_results_df(efit_list[[i]])
  }),
  paste0("efit_", names(efit_list), "_results_df") # Informative names
)

top_DE_entrezIDs <- function(df, direction = "up", min_genes = 5) {
  msg <- NULL
  
  # Initial filter
  if (direction == "up") {
    filtered_genes <- df %>% filter(logFC > 0.58 & adj.P.value < 0.05)
  } else if (direction == "down") {
    filtered_genes <- df %>% filter(logFC < -0.58 & adj.P.value < 0.05)
  } else {
    stop("Invalid direction. Choose 'up' or 'down'")
  }
  
  # If too few genes, relax criteria
  if (nrow(filtered_genes) < min_genes) {
    msg <- glue::glue(
      "⚠️ {nrow(filtered_genes)} genes passed conservative DE filter for direction '{direction}'. ",
      "Relaxed filter used (logFC > 0 & adj.P < 0.1)."
    )
    
    if (direction == "up") {
      filtered_genes <- df %>% filter(logFC > 0 & adj.P.value < 0.1)
    } else {
      filtered_genes <- df %>% filter(logFC < 0 & adj.P.value < 0.1)
    }
  }
  
  # Still empty?
  if (nrow(filtered_genes) == 0) {
    msg <- glue::glue(
      "❌ No DE genes found for direction '{direction}', even after relaxing thresholds."
    )
    return(list(entrez_ids = character(0), message = msg))
  }
  
  # Try mapping to Entrez IDs
  entrez_ids <- tryCatch({
    mapIds(
      annotation_obj,
      keys = filtered_genes$ensembleID,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  }, error = function(e) {
    msg <<- glue::glue(
      "❌ No DE genes found mapped to Entrez IDs for direction '{direction}': {conditionMessage(e)}"
    )
    return(character(0))
  })
  
  return(list(entrez_ids = na.omit(entrez_ids), message = msg))
}

# Apply function dynamically for each contrast and store results
top_DE_entrezIDs_list <- setNames(
  lapply(efit_results_list, function(df) {
    list(
      up = top_DE_entrezIDs(df, "up"),
      down = top_DE_entrezIDs(df, "down")
    )
  }),
  paste0("top_DE_entrezIDs_", names(efit_results_list))
)


universe_entrez <- mapIds(
  annotation_obj,
  keys = rownames(DGE.NOIseqfilt),  # Ensembl universe
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Rename some objs to work with exisiting .Rmd code
dge_list_edgeR <- DGE.edgeRfilt
dge_list_raw <- DGE
dge_list_NOISeq <- DGE.NOIseqfilt
efit_list <- efit_list
efit_results_dfs <- efit_results_list
entrez_ids_list <- top_DE_entrezIDs_list


