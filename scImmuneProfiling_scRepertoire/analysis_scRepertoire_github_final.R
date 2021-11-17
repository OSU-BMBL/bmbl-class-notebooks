rm( list=ls() )
graphics.off()
# library(devtools)
# devtools::install_github("ncborcherding/scRepertoire")
library(scRepertoire)
library(Seurat)
library(circlize)
library(scales)


# ------------------- Loading and Processing Data ------------------- 

# Loading and Processing Data
# Borcherding et al., (2021), Communications Biology, 4: 122

# Load Seurat Object of scRNA-seq Data

seurat <- get(load("seurat2.RData"))

# Load scTCR-seq Data

# S1 <- read.csv(".../Sample1/outs/filtered_contig_annotations.csv")
# S2 <- read.csv(".../Sample2/outs/filtered_contig_annotations.csv")
# S3 <- read.csv(".../Sample3/outs/filtered_contig_annotations.csv")
# S4 <- read.csv(".../Sample4/outs/filtered_contig_annotations.csv")
# S5 <- read.csv(".../Sample5/outs/filtered_contig_annotations.csv")
# S6 <- read.csv(".../Sample6/outs/filtered_contig_annotations.csv")
# 
# contig_list <- list(S1, S2, S3, S4, S5, S6)

data("contig_list") #the data built into scRepertoire

head(contig_list[[1]])

# Combining and Annotating the Data

combined <- combineTCR(contig_list, 
  samples = c("PY", "PY", "PX", "PX", "PZ","PZ"), 
  ID = c("P", "T", "P", "T", "P", "T"), cells ="T-AB")

# Adding Additional Variables

example <- addVariable(combined, name = "batch", 
  variables = c("b1", "b1", "b2", "b2", "b2", "b2"))

names(example)
head(example[[1]])

# Subsetting

subset <- subsetContig(combined, name = "sample", 
  variables = c("PX", "PY"))

names(subset)
head(subset[[1]])


# ------------------- Clonotype Analysis ------------------- 

# Quantify Clonotypes

quantContig(combined, cloneCall="gene+nt", scale = TRUE)

quantContig(combined, cloneCall="gene+nt", chain = "TRA")

quantContig(combined, cloneCall="gene+nt", chain = "TRB")

quantContig(combined, cloneCall="gene+nt", chain = "TRA", scale = TRUE)

quantContig(combined, cloneCall="gene+nt", chain = "both", scale = TRUE)

quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
  scale = TRUE, exportTable = TRUE)
quantContig_output

# Track Clonotypes

compareClonotypes(combined, numbers = 10, samples = c("PX_P", "PX_T"), 
  cloneCall="aa", graph = "alluvial")

# Visualize Gene Usage

vizGenes(combined, gene = "V", chain = "TRB", plot = "bar", 
  order = "variance", scale = TRUE)

# Clonal Space Homeostasis

clonalHomeostasis(combined, cloneCall = "gene")

clonalHomeostasis(combined, cloneCall = "aa")

# Clonal Proportion

clonalProportion(combined, cloneCall = "gene") 

clonalProportion(combined, cloneCall = "nt") 

# Overlap Analysis

clonalOverlap(combined, cloneCall = "gene+nt", method = "morisita")

# Diversity Analysis

clonalDiversity(combined, cloneCall = "gene", groupBy = "sample", x.axis = "ID")

# Scatter Compare

scatterClonotype(combined, cloneCall ="gene", 
  x.axis = "PY_P", y.axis = "PY_T", dot.size = "total", graph = "proportion")


# ------------------- Interacting with Seurat ------------------- 

# UMAP Visualization

DimPlot(seurat, label = TRUE) + NoLegend()
table(Idents(seurat))

DimPlot(seurat, group.by = "Type")

# Adding Clonotype Information

seurat <- combineExpression(combined, seurat, cloneCall="gene+nt")
seurat@meta.data[1:5,]

# UMAP Visualization by Clonotype

DimPlot(seurat, group.by = "cloneType")

seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CAVNGGSQGNLIF_CSAEREDTDTQYF", "NA_CATSATLRVVAEKLFF"))
DimPlot(seurat, group.by = "highlight")

# Clonotype vs. Cell Clusters

occupiedscRepertoire(seurat, x.axis = "cluster")

# Clonotype Tracking: Alluvial Plot

alluvialClonotypes(seurat, cloneCall = "gene", 
  y.axes = c("Patient", "cluster", "Type"), 
  color = "cluster") 

# Clonotype Tracking: Chord Diagram 

circles <- getCirclize(seurat, groupBy = "cluster")

grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)

circlize::chordDiagram(circles, self.link = 1, grid.col = grid.cols)
