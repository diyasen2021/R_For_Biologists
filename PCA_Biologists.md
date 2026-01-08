# PCA for Biologists (RNA-seq context, R-friendly)

## What problem PCA solves

In RNA-seq experiments, each sample has expression values for **thousands of genes**. This makes it difficult to directly see patterns in the samples such as:

- Do biological replicates cluster together?
- Do treated and control samples separate?
- Are there batch effects or outlier samples?

**Principal Component Analysis (PCA)** simplifies this complexity by reducing thousands of gene expression measurements into just a few **summary axes** also called **principal components** that capture the main differences between samples.

---

## What PCA does

PCA takes the following as input
- **Samples** → the things you want to compare  
- **Genes** → features used to compare samples  

- Each sample is represented as a long list of gene expression values.
- PCA looks at **all genes together** to determine how samples differ overall.
- PCA compresses the gene matrix into a single axis (principal component) by combining the genes into weighted summaries 
- Each principal component is therefore not a single gene, but a weighted combination of many genes, with genes that show large, consistent differences between samples contributing more strongly.
- PCA ranks these components by how much variation they capture, keeping only the most informative ones (for example, the first 5 PCs)
  

---

## What each point in a PCA plot means

In a PCA plot:

- **Each point = one sample**
- Samples close together have **similar gene expression profiles**
- Samples far apart differ in **many genes**

The distance between points reflects how different the samples are **across all genes**, not just one or two.

What PC1 and PC2 mean

### PC1 (Principal Component 1)
The direction that explains the **largest source of variation** between samples.

### PC2 (Principal Component 2)
The second largest source of variation, independent of PC1.

The percentage shown on each axis indicates how much of the **total variation** is explained by that component.

**Example:**

- `PC1 (35%)` means that 35% of all sample-to-sample variation is captured along that axis.

---

## How samples are placed on the plot

PCA:

1. Identifies genes that vary the most across samples  
2. Uses these genes to define new axes (principal components)  
3. Calculates a **score** for each sample along each axis  
4. Plots these scores as coordinates  

The plot is **entirely data-driven**. PCA does **not** know which samples are controls or treatments.

---

## Why PCA is useful in RNA-seq analysis

PCA is mainly used for **exploratory analysis and quality control**:

- Check if biological replicates cluster together
- See whether samples separate by condition
- Detect batch effects
- Identify outlier samples

If samples separate by condition on PC1 or PC2, it suggests that the condition has a strong effect on gene expression.

---

## Important note about preprocessing (very important in R)

PCA should **not** be run on raw counts.

In R, PCA is typically performed on:

- **Normalized counts**
- **Log-transformed or variance-stabilized data**  
  (e.g. `vst()` or `rlog()` from **DESeq2**)

This prevents highly expressed genes or sequencing depth from dominating the PCA.

---

## How PCA is typically done in R (conceptual)

Although R handles the mathematics internally, conceptually PCA in R:

- Takes a **samples × genes** matrix
- Finds major axes of variation
- Returns coordinates (scores) for each sample
- Plots PC1 vs PC2

You interpret the plot **biologically**, not statistically.

---

## Minimal PCA example in R using DESeq2

This is a **minimal, clean example**, showing only the essential steps.


### 1. Load libraries

```r
library(DESeq2)
library(ggplot2)
```

### 2. Download some data
You can download the processed count and metadata files directly from GitHub:

Counts data (genes × samples)
https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_counts_cerebellum.csv
 
Sample metadata (conditions, etc.)
https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_coldata_cerebellum.cs

### 3. Import data into R

```
counts <- read.csv("GSE96870_counts_cerebellum.csv", row.names=1)
coldata <- read.csv("GSE96870_coldata_cerebellum.csv", row.names=1)
```

### 4. Create a DESeq2 object

```r
dds <- DESeqDataSetFromMatrix(
  countData = counts,      # genes x samples matrix
  colData   = coldata,     # sample information
  design    = ~ condition  # experimental condition
)
```

### 5. Run DESeq2 preprocessing

```r
dds <- DESeq(dds)
```

### 6. Apply variance-stabilizing transformation

```r
vsd <- vst(dds, blind = TRUE)
```

### 7. Perform PCA and plot

```r
plotPCA(vsd, intgroup = "condition")
```

---

## 8. How to interpret this PCA plot

- **Each point = one sample**
- Samples close together → similar gene expression profiles
- Separation along PC1 or PC2 → major transcriptomic differences
- Coloring by `condition` helps assess whether biology explains the variation

Axis labels typically look like:

```
PC1: 38% variance
PC2: 21% variance
```

---

## Important Notes

- **PCA is unsupervised**
  - `condition` is not used to calculate PCs
  - It is only used for coloring the plot

- PCA is for:
  - Quality control
  - Detecting batch effects
  - Checking replicate consistency

- PCA is **not** used to identify differentially expressed genes

---

## One-sentence summary

> PCA reduces thousands of gene expression values into a few components that summarize the main differences between RNA-seq samples, allowing us to visualize similarity, detect batch effects, and assess experimental quality.
