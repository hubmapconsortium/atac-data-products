#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(ArchR)
    library(Seurat)
    library(ShinyCell2)
})

############################################
## 1. Parse command-line options
############################################

option_list <- list(
    make_option(c("--input"), type="character", help="Path to input fragments.tsv.gz"),
    make_option(c("--outdir"), type="character", help="Output directory"),
    make_option(c("--threads"), type="integer", default=4, help="Threads [default=4]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) | is.null(opt$outdir)) {
    stop("\n❌ ERROR: Missing required --input or --outdir argument\n")
}

input_file <- normalizePath(opt$input)
output_dir <- normalizePath(opt$outdir, mustWork = FALSE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message(">>> INPUT FILE: ", input_file)
message(">>> OUTPUT DIR: ", output_dir)
message(">>> THREADS: ", opt$threads)

############################################
## 2. ArchR Setup
############################################

addArchRGenome("hg38")
addArchRThreads(threads = opt$threads)

sample_name <- tools::file_path_sans_ext(basename(input_file))
inputFiles <- c(input_file)
names(inputFiles) <- sample_name

############################################
## 3. Create Arrow Files
############################################

message(">>> Creating Arrow files ...")

ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    minTSS = 1,
    minFrags = 1,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)

############################################
## 4. Create ArchR Project
############################################

proj_dir <- file.path(output_dir, "ArchR_Project")
message(">>> Creating ArchRProject at: ", proj_dir)

archr_proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = proj_dir,
    copyArrows = TRUE
)

if (!inherits(archr_proj, "ArchRProject")) {
    stop("❌ ERROR: Failed to create ArchR project (ArchRProj is NULL or invalid)")
}

############################################
## 5. LSI / Clustering
############################################

archr_proj <- addIterativeLSI(
    archr_proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(resolution = c(0.2), n.start = 5),
    varFeatures = 25000,
    dimsToUse = 1:5
)

if (!inherits(archr_proj, "ArchRProject")) {
    stop("❌ ERROR: IterativeLSI returned an invalid ArchR project")
}

archr_proj <- addClusters(
    archr_proj, 
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters"
)

## Overwrite clusters with two pseudo groups
n <- nCells(archr_proj)
archr_proj$Clusters <- factor(c(rep("C1", ceiling(n/2)), rep("C2", floor(n/2))))

############################################
## 6. UMAP
############################################

archr_proj <- addUMAP(
    archr_proj,
    reducedDims = "IterativeLSI",
    nNeighbors = 2,
    force = TRUE
)

## Fix cellColData row mismatch
valid_cells <- intersect(archr_proj$cellNames, rownames(archr_proj@embeddings$UMAP$df))
archr_proj@cellColData <- archr_proj@cellColData[valid_cells, , drop = FALSE]

############################################
## 7. Save project
############################################

archr_proj <- saveArchRProject(
    archr_proj,
    outputDirectory = proj_dir,
    load = TRUE
)

############################################
## 8. Pseudo-peak generation
############################################

featureDF <- getFeatures(archr_proj, useMatrix = "TileMatrix")
myPeaks <- GenomicRanges::GRanges(
    seqnames = featureDF@seqnames,
    ranges   = IRanges::IRanges(
        start = featureDF@ranges@start,
        end   = featureDF@ranges@start + 499
    ),
    name = featureDF$idx
)

archr_proj <- addPeakSet(archr_proj, peakSet = myPeaks)
archr_proj <- addPeakMatrix(archr_proj)

############################################
## 9. ShinyCell2 Setup
############################################

message(">>> Preparing ShinyCell2 metadata...")

scConf2 <- createConfig(archr_proj)

## Remove unwanted QC columns
to_del <- c("ReadsInTSS", "ReadsInPromoter", "ReadsInBlacklist",
            "NucleosomeRatio", "nMultiFrags", "nMonoFrags", "nFrags",
            "nDiFrags", "ReadsInPeaks")

to_del_existing <- to_del[to_del %in% colnames(scConf2$meta)]
scConf2 <- delMeta(scConf2, to_del_existing)

############################################
## 10. Generate ShinyCell2 Output
############################################

shiny_dir <- file.path(output_dir, "ShinyCell2")
dir.create(shiny_dir, showWarnings = FALSE, recursive = TRUE)

message(">>> Generating ShinyCell2 files at: ", shiny_dir)

makeShinyFiles(
    archr_proj,
    scConf = scConf2,
    dimred.to.use = "UMAP",
    bigWigGroup = c("Clusters"),
    shiny.prefix = "scATAC",
    shiny.dir = shiny_dir,
    default.gene1 = "IRF1",
    default.multigene = NA,
    default.dimred = "UMAP",
    chunkSize = 500
)

makeShinyCodes(
    shiny.title = "DEMO of Heart Data: scATAC shinycell2",
    shiny.prefix = c("scATAC"),
    shiny.headers = c("scATAC-seq Data"),
    shiny.dir = shiny_dir
)

message("\n============================================")
message("   ✅ FINISHED! ShinyCell2 Output Saved To:")
message("   ", shiny_dir)
message("============================================\n")
