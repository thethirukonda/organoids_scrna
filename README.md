# Single-Cell RNA Sequencing Analysis of Human Cerebral Organoids

Comprehensive scRNA-seq analysis workflow for characterizing neural development in human cerebral organoids, integrating two independent datasets to identify cell types, developmental trajectories, and gene regulatory mechanisms underlying forebrain development.

This tutorial is based on -  [Zhisong He's tutorial for scrna](https://github.com/quadbio/scRNAseq_analysis_vignette/blob/master/README.md)

i have made some changes - especially using a intergrated dataset where i combine2 different datasets of human brain organoids


## Overview

This analysis integrates single-cell transcriptomic data from two  human cerebral organoid datasets to comprehensively characterize neural development. The workflow combines multiple validation approaches including cluster annotation, pseudotemporal trajectory analysis, transcriptome similarity comparisons, reference-based label transfer, and pathway enrichment analysis to build a robust map of forebrain development from neural progenitors to mature neurons.

## Project Goals

- **Batch correction & integration:** Remove technical variation between DS1 and DS2 while preserving biological signal using Harmony, enabling joint analysis of ~20,000 cells
- **Cell type discovery & annotation:** Identify and characterize biologically distinct cell populations through unsupervised clustering and marker gene analysis, mapping the full spectrum of neural development from progenitors to mature neurons
- **Dorsal-ventral specification:** Distinguish fundamental forebrain developmental axis (dorsal excitatory vs ventral GABAergic lineages) through regional clustering and lineage-specific differential gene expression analysis
- **Developmental pseudotemporal ordering:** Map dorsal developmental progression (progenitor → intermediate progenitor → mature neuron) using diffusion pseudotime, revealing temporal gene expression dynamics across developmental stages
- **Gene expression characterization:** Identify differentially expressed genes and stage-specific transcriptional programs driving progenitor proliferation, intermediate progenitor transition, and neuronal maturation
- **Multi-method annotation validation:** Confirm cell type assignments through four independent approaches (transcriptome similarity, cell-level correlation, Seurat label transfer, marker gene consistency) achieving >90% cross-method agreement
- **Functional interpretation:** Use gene set enrichment analysis to connect differential gene expression to biological pathways (cell cycle, synaptogenesis, excitatory/inhibitory neurotransmission) revealing molecular mechanisms underlying development
- **Data accessibility:** Export annotated data (H5AD format) for downstream Python-based analyses (e.g., RNA velocity with scVelo)

## Data Description

### Input Datasets

**Dataset 1 (DS1):** Primary cerebral organoid dataset
- Genes: ~20,000
- Cells: ~11,000 (pre-filtering)
- Platform: 10x Genomics 3' scRNA-seq

**Dataset 2 (DS2):** Secondary cerebral organoid dataset  
- Genes: ~20,000
- Cells: ~12,000 (pre-filtering)
- Platform: 10x Genomics 3' scRNA-seq

### Quality Control Filtering

Applied consistent filtering to both datasets:
- **Min genes per cell:** 200 genes
- **Max genes per cell:** 6,000 genes
- **Mitochondrial content threshold:** <5%
- **Min cells expressing gene:** 3

**Final cell count:** ~20,000 high-quality cells (combined after QC)

## Analysis Workflow

### 1. Data Loading and Initial Quality Control

Read 10x Count matrices for both datasets:

```R
ds1_counts <- Read10X(data.dir = "data/DS1")
ds2_counts <- Read10X(data.dir = "data/DS2")
```

Created Seurat objects with initial quality thresholds:
- Minimum 3 cells expressing each gene
- Minimum 200 genes per cell

Calculated QC metrics:
- **nFeature_RNA:** Genes detected per cell
- **nCount_RNA:** UMI counts per cell  
- **percent.mt:** Mitochondrial gene percentage

**Rationale:** Removes low-quality libraries, ambient contamination, and identifies stressed/dying cells.

### 2. Quality Filtering

Applied consistent filtering criteria to both datasets:

```R
subset(nFeature_RNA > 200 & nFeature_RNA < 6,000 & percent.mt < 5)
```

Removes:
- Doublets (exceptionally high UMI/gene counts)
- Empty droplets (very low counts)
- Low-quality cells (few genes detected)
- Stressed/dying cells (high mitochondrial %)

### 3. Normalization and Feature Selection

Applied independent normalization to each dataset:

**Per-dataset processing:**
```R
NormalizeData()        # Log-normalization (CPM scaling)
FindVariableFeatures() # 2,000 highly variable genes (HVGs)
ScaleData()            # Z-score normalization across genes
```

**Rationale for 2,000 HVGs:** Balance between capturing biological variation and computational efficiency. Tutorial used 3,000 but 2,000 is more standard.

**For merged dataset:** 3,000 HVGs (accounts for increased variation from combining two sources)

### 4. Visualization of Normalization

Generated variable feature plots showing:
- Relationship between mean expression and variance
- Top 20 most variable genes highlighted
- Confirms removal of mitochondrial and ribosomal genes

### 5. Batch Effect Visualization (Pre-Correction)

Merged datasets without integration to visualize batch effects:

```R
merged <- merge(ds1, ds2) %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50) %>%
    RunUMAP(dims = 1:20)
```

**Observation:** Clear separation by dataset (`orig.ident`), confirming batch effect requiring correction.

**Marker genes checked:** FOXG1 (telencephalon), EMX1 (dorsal), DLX2 (ventral), LHX9 (non-telencephalon)

### 6. Batch Correction with Harmony

Applied Harmony batch correction to remove dataset-specific effects:

```R
RunHarmony(merged, group.by.vars = "orig.ident", 
           dims.use = 1:20, max.iter.harmony = 50)
```

**Parameters:**
- Batch variable: Dataset source (DS1 vs DS2)
- Iterations: 50 (default convergence)
- Dimensions: First 20 PCs

**Why Harmony?** Fast, scalable, preserves biological variation. Alternative methods (Seurat Integration, Combat) exist but Harmony provides good balance.

**Validation:** Post-correction UMAP shows mixed colors across space, confirming successful batch removal.

### 7. Dimensionality Reduction

Applied to harmony-corrected data:

**PCA:** 50 principal components
- Identifies major sources of variance
- Used for downstream dimension reduction

**Elbow Plot Analysis:** 
- Determined when variance explained plateaus
- Chose 20 PCs for downstream analysis (balances detail vs noise)

**UMAP:** Using harmony reduction, dimensions 1:20
- 2D visualization for cluster visualization
- Non-linear dimension reduction preserves local structure

**t-SNE:** Alternative 2D visualization
- Compared with UMAP for consistency
- Both show similar cluster organization

**PCA Heatmaps:** 
- Visualized top genes driving first 9-18 PCs
- Confirmed expected developmental genes driving variation

### 8. Clustering

**Method:** Louvain algorithm on k-nearest neighbor graph

```R
FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
FindClusters(resolution = 0.6)
```

**Two-step resolution exploration:**
1. Resolution 0.6: Balanced granularity
2. Resolution 1.0: Finer-grained clusters for comparison

**Result:** 14-16 clusters, representing distinct transcriptomic states

### 9. Cell Type Annotation Strategy

#### Marker Gene Selection

Used comprehensive marker gene panel:

**Telencephalic Identity Markers:**
- FOXG1 (telencephalic specification)
- EMX1 (dorsal telencephalon)
- DLX2, DLX5 (ventral telencephalon)
- NKX2-1 (ventral forebrain, GABAergic)
- PAX6 (progenitor, dorsal cortex)

**Developmental Stage Markers:**
- MKI67 (proliferation, G2/M phase)
- NES (intermediate filament, neural progenitors)
- DCX (neuroblast, immature neurons)
- SOX2 (radial glia, stem cells)
- EOMES (intermediate progenitors, TBR2)
- NEUROD6 (mature excitatory neurons, NeuroD)

**Non-Telencephalic Markers (Contamination Detection):**
- OTX2, LHX9, TFAP2A (midbrain/hindbrain)
- RELN, HOXB family (non-cortical)
- RSPO3 (mesenchymal contamination)
- GATA3, AIF1, VIT, PTPRC (blood/immune cells)

**Note:** Markers derived from literature; for your own datasets, perform literature review or use PanglaoDB/CellMarker.

#### Annotated Cell Types

| Cluster | Cell Type | Key Features |
|---------|-----------|--------------|
| 0 | Mature dorsal excitatory neurons | DCX+, NEUROD6+ neurons |
| 1 | Cycling progenitors | MKI67+ highly proliferative |
| 2 | Ventral interneurons | DLX2+ developing inhibitory |
| 3 | Ventral progenitors | DLX5+, NKX2-1+ expanding pool |
| 4 | GABAergic interneurons | GABAergic neurotransmitter phenotype |
| 5 | Hippocampal/cortical hem | Regional specification markers |
| 6 | G2/M cycling cells | Active S/G2 cell cycle phase |
| 7 | Intermediate progenitors | EOMES+, developmental transition |
| 8 | SST+ interneurons | Somatostatin-expressing subset |
| 9 | Mixed/transitional cells | Between developmental states |
| 10 | Dorsal progenitors | PAX6+, EMX1+, high proliferation |
| 11 | S-phase cycling cells | S phase cell cycle phase |
| 12 | Choroid plexus/roof plate | Regional specialization |
| 13 | G1/S cycling progenitors | G1/S cell cycle phase |
| 14+ | Contaminating cells | Blood, mesenchymal (removed) |

**Contamination Removal:** Clusters 14-15 showing blood markers removed from downstream analysis.

### 10. Validation Approaches

Multiple independent methods confirm cell type assignments:

#### Method 1: Transcriptome Similarity (Cluster Level)

Calculated Spearman correlation between average cluster expression vs reference cell types:

```R
corr2ref_cl <- cor(avg_expr_yours, avg_expr_ref, method = "spearman")
```

**Procedure:**
- Computed average gene expression per cluster
- Computed average gene expression per reference cell type
- Calculated correlation across ~15,000 common genes
- Clustered heatmap reveals which of your clusters match reference types

**Results:** Correlations 0.94-0.98, high diagonal, confirming annotations

**Interpretation:** Strong correlation validates cluster identities match reference cell types

#### Method 2: Transcriptome Similarity (Cell Level)

Assigned each individual cell to best-matching reference cell type:

```R
ranked_expr_ref <- apply(avg_expr_ref, 2, rank)
ranked_expr_yours <- rank_matrix(integrated_clean@assays$RNA@data)
corr2ref_cell <- corSparse(ranked_expr_yours, ranked_expr_ref)
```

**Advantage:** Finer resolution than cluster-level
- Identifies transitional cells
- Reveals within-cluster heterogeneity

**Finding:** While cluster-level correlations high, individual cells exist along continuum rather than discrete states

#### Method 3: Seurat Label Transfer

Reference-based label transfer learning:

```R
anchors <- FindTransferAnchors(reference = seurat_ref, query = integrated_clean)
predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$celltype)
```

**Process:**
- Find anchor cells between query and reference
- Transfer pre-computed reference labels
- Compare agreement with manual annotation

**Results:**
- Homogeneous clusters: >70% assignment to single cell type
- Heterogeneous clusters: Mixed predictions indicate transitional populations
- Overall agreement with manual annotations: >90%

**Interpretation:** Transitional cells represent differentiation in progress

#### Method 4: Marker Gene Consistency

Final validation: examined comprehensive marker genes

```R
cl_markers <- FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
```

- Identified top 10 significant markers per cluster
- Verified cell-type-specific expression patterns
- Cross-checked against published organoid literature

**Result:** All major populations showed expected marker genes

### 11. Pseudotemporal Trajectory Analysis (Dorsal Lineage)

Focused on dorsal forebrain development: progenitors → intermediate progenitors → mature neurons

#### Subset Dorsal Clusters

Selected developmental trajectory clusters:
- Cluster 0: Mature excitatory neurons (endpoint)
- Cluster 7: Intermediate progenitors (transition)
- Cluster 10: Dorsal progenitors (starting point)

#### Cell Cycle Regression for Trajectory

**Critical step:** Removed cell cycle genes from variable features

```R
VariableFeatures(integrated_dorsal) <- setdiff(VariableFeatures(integrated_dorsal), 
                                                unlist(cc.genes))
```

**Rationale:** Cell cycle dominates transcriptional variance in proliferating progenitors. Even though we don't regress out values, removing from PCA reveals underlying developmental signal.

**Result:** Clusters less mixed, trajectories clearer after cell cycle removal

#### Diffusion Pseudotime (DPT) Computation

Computed pseudotime using diffusion map distance:

```R
dm <- DiffusionMap(Embeddings(integrated_dorsal, "pca")[, 1:20])
dpt <- DPT(dm)
integrated_dorsal$dpt <- rank(dpt$dpt)
```

**Interpretation:**
- Pseudotime = computational measure of developmental progression
- Cells ordered by diffusion distance through transcriptional space
- Ranking converts to 1-N pseudotime values

#### Developmental Gene Dynamics

Plotted gene expression across pseudotime with LOESS smoothing:

**GLI3 (Early Progenitor Marker):**
- Highest at trajectory start
- Rapidly decreases with increasing pseudotime
- Confirms it marks "starting state"

**EOMES (Intermediate Progenitor Marker):**
- Shows "pulse" or wave pattern
- Peaks mid-trajectory (~rank 750)
- Drops toward trajectory end
- Identifies intermediate progenitor stage

**NEUROD6 (Mature Neuron Marker):**
- Increases throughout trajectory
- Stays high at trajectory end
- Confirms neurons as final differentiation state

**Sparsity/Dropout:** Thick band of zeros at Y=0
- Technical artifact of scRNA-seq, not biological
- Many cells don't express gene (mRNA capture limitation)
- Legitimate biological absence can't be distinguished

**LOESS smoothing:** Blue line reveals biological trend through noisy data points

### 12. Differential Gene Expression Analysis

Identified genes driving developmental transitions and lineage specification.

#### A. Progenitors vs Mature Neurons

Most dramatic developmental transition:

```R
deg_prog_vs_neuron <- FindMarkers(
    ident.1 = "10",    # Dorsal progenitors
    ident.2 = "0",     # Mature neurons
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = FALSE
)
```

**Statistics:**
- Total DEGs (padj < 0.05): ~2,500 genes
- Upregulated in progenitors: ~1,000 genes
- Upregulated in neurons: ~1,500 genes

**Top Progenitor Genes:** MKI67, PCNA, TOP2A
- Cell division and proliferation machinery
- Cell cycle checkpoints and regulation

**Top Neuron Genes:** SYN1, SNAP25, VAMP2
- Synaptic vesicle transport
- Synaptic transmission machinery
- Neurotransmitter release

**Visualization:**
- Volcano plot: significance vs effect size
- Heatmap of top genes shows clean separation between clusters

#### B. Dorsal vs Ventral Lineages

Fundamental developmental axis (forebrain organization):

**Region Definition:**
- Dorsal: Clusters 0, 7, 10, 1, 11 (excitatory neurons, dorsal progenitors)
- Ventral: Clusters 2, 4, 8 (GABAergic neurons, ventral progenitors)

```R
deg_dorsal_vs_ventral <- FindMarkers(
    ident.1 = "Dorsal",
    ident.2 = "Ventral",
    test.use = "wilcox"
)
```

**Statistics:**
- Total DEGs: ~2,000 genes
- Dorsal markers: ~1,200 (glutamatergic signaling)
- Ventral markers: ~800 (GABAergic signaling)

**Top Dorsal Genes:** NEUROD1, EMX1, CUX1, CUX2
- Excitatory neurotransmission
- Cortical lamination and patterning
- Forebrain-specific transcription factors

**Top Ventral Genes:** GAD1, GAD2, DLX2, NKX2-1
- GABAergic neurotransmitter synthesis
- Interneuron specification
- Ventral forebrain development

**Spatial Observation:** Clean UMAP boundaries separate dorsal and ventral regions, confirming correct regional assignments

**Biological Significance:** 
- Dorsal lineage → excitatory neurons (cortex, thinking/planning)
- Ventral lineage → inhibitory neurons (regulation)
- This axis is fundamental feature of forebrain development

#### C. Three-Way Maturation Trajectory

Stage-specific expression across Progenitor → IP → Neuron progression:

**Stage Definition:**
- Progenitor: Clusters 9, 10, 11 (cycling, high proliferation)
- IP (Intermediate Progenitor): Cluster 7 (transitional stage)
- Neuron: Cluster 0 (mature, post-mitotic)

**Three Pairwise Comparisons:**

1. **Progenitor vs IP:**
   - Progenitor-high genes: SOX2, PCNA (stemness, proliferation)
   - IP-high genes: NEUROG1, NEUROG3 (neurogenic transcription factors)

2. **IP vs Neuron:**
   - IP-high genes: EOMES, NEUROD2 (developmental transition)
   - Neuron-high genes: SYN1, SNAP25 (synaptic maturation)

3. **Progenitor vs Neuron:**
   - Most dramatic changes
   - Combines effects of both earlier transitions

**Stage-Specific Genes:**
- **IP-specific** (high in IP only): ~1,951 genes
  - Unique to IP developmental switch
  - Identify IP-specific transcriptional programs
- **Progressive** (change across all transitions): ~1,200 genes
  - Continuously up/down-regulated
  - Drive differentiation across all stages

### 13. Gene Set Enrichment Analysis (GSEA)

Interpreted differential expression in context of biological pathways and processes.

**Method:** Gene Ontology (GO) enrichment using `enrichGO` from clusterProfiler

#### GSEA 1: Progenitor vs Neuron Pathways

**Progenitor-Enriched Pathways:**
- Cell cycle regulation
- DNA replication and repair
- Mitochondrial biogenesis
- Metabolic biosynthetic pathways
- Organelle assembly

**Interpretation:** Progenitors have high metabolic demand and continuous division. Must replicate DNA, produce organelles, maintain high energy.

**Neuron-Enriched Pathways:**
- Synaptic transmission
- Neurotransmitter signaling
- Calcium signaling (synaptic plasticity)
- Axon guidance and morphogenesis
- Long-term potentiation (learning/memory)
- Myelination and myelin formation

**Interpretation:** Neurons establish functional circuits. Synaptic connections enable information transfer. Plasticity allows learning.

#### GSEA 2: Dorsal vs Ventral Pathways

**Dorsal-Enriched Pathways:**
- Forebrain development (TelencephalicDevelopment)
- Cortical patterning and regional identity
- Glutamate signaling and excitatory transmission
- Dendritic morphogenesis (building neural structure)
- Neuronal projection morphogenesis

**Interpretation:** Dorsal forebrain specialized for glutamatergic excitation, building cortical circuits for cognition.

**Ventral-Enriched Pathways:**
- GABAergic neuron development
- Interneuron specification and migration
- GABAergic synapse organization
- Synaptic inhibition

**Interpretation:** Ventral forebrain specialized for GABAergic inhibition, generating interneuron diversity.

#### GSEA 3: IP-Specific Pathways

Identified pathways unique to intermediate progenitor stage:

Genes enriched in IP (but not progenitors or neurons):

**IP-Specific Enriched Pathways:**
- Neural differentiation and cell fate commitment
- Developmental signaling (BMP, Wnt, Notch)
- Transcriptional regulation and chromatin remodeling
- Neuroblast differentiation

**Interpretation:** IP stage represents critical developmental switch-point with distinct molecular programs preparing transition to neuron.

#### GSEA Visualization

For each analysis:
- **Dot plots:** Top 20 enriched pathways
  - Size = gene count in pathway
  - Color = p-value (red = significant)
  - Shows pathway enrichment landscape
- **Bar plots:** Top 10 pathways for clarity
- **CSV exports:** Full results with p-values, gene counts

**Combined Visualization:** Multi-panel plot comparing:
- Progenitor vs Neuron pathways (2 panels)
- Dorsal vs Ventral pathways (2 panels)
- IP-Specific pathways (1 panel)

### 14. Data Export for Downstream Analysis

Converted final Seurat object to H5AD format for Python-based analysis:

```R
adata <- AnnData(
    X = t(as.matrix(GetAssayData(integrated_clean, layer = "data"))),
    obs = integrated_clean@meta.data,
    var = data.frame(row.names = rownames(integrated_clean))
)

# Add all embeddings
adata$obsm <- list(
    X_umap = Embeddings(integrated_clean, "umap"),
    X_harmony = Embeddings(integrated_clean, "harmony"),
    X_pca = Embeddings(integrated_clean, "pca")
)

adata$write_h5ad("integrated_clean.h5ad")
```

**Output:** `integrated_clean.h5ad`
- Contains: Expression data, metadata, all embeddings
- Format: AnnData (scanpy-compatible)
- Use case: RNA velocity analysis with scVelo in Python

## Software and Computational Requirements

### R Packages Used

```R
library(dplyr)          # Data wrangling
library(Seurat)         # Single-cell analysis (v5.0+)
library(patchwork)      # Figure composition
library(ggplot2)        # Publication-quality plotting
library(ggrepel)        # Smart text labels for plots
library(pheatmap)       # Heatmap visualization
library(harmony)        # Batch effect correction
library(destiny)        # Diffusion maps for pseudotime
library(org.Hs.eg.db)   # Human gene ID mappings
library(clusterProfiler)  # Gene set enrichment
library(presto)         # Fast Wilcoxon testing
library(qlcMatrix)      # Sparse matrix correlations
library(anndata)        # H5AD format I/O
```


### Parameter Settings

| Step | Parameter | Setting | Purpose |
|------|-----------|---------|---------|
| Normalization | Method | Log-normalization | Stabilize variance |
| HVGs (individual) | nfeatures | 2,000 | Balance variance capture |
| HVGs (merged) | nfeatures | 3,000 | Account for dataset variation |
| PCA | npcs | 50 | Capture major variance |
| Harmony | dims.use | 1:20 | Remove technical noise |
| Harmony | max.iter | 50 | Default convergence |
| UMAP/neighbors | dims | 1:20 | Balances detail vs speed |
| Clustering | resolution | 0.6-1.0 | Explore granularity |
| FindMarkers | min.pct | 0.25 | Gene expressed in >25% |
| FindMarkers | logfc.threshold | 0.25 | Minimum effect size |
| Pseudotime | method | DPT | Robust trajectory inference |
| GSEA | pvalueCutoff | 0.05 | Standard significance |
| GSEA | ont | BP | Biological Processes GO |



## File Structure and Outputs

# R mark down 

scrna_brain

### Data Files
```
data/
├── DS1/                    # Dataset 1 10x output
│   ├── matrix.mtx
│   ├── barcodes.tsv
│   ├── features.tsv
│   └── ref_seurat_obj.rds  # Reference for validation
└── DS2/                    # Dataset 2 10x output
    ├── matrix.mtx
    ├── barcodes.tsv
    └── features.tsv
```

### Results - Figures (PDF)
```
All figures are saved in `plots_and_graphs/`

### **Quality Control & Preprocessing**

| File | Description |
|------|-------------|
| `ds1violinplot.png` | QC metrics for Dataset 1 (nFeature, nCount, %MT) |
| `ds2violinplot.png` | QC metrics for Dataset 2 |
| `ds1scatterplot.png` | Feature-count relationship DS1 |
| `ds2scatterplot.png` | Feature-count relationship DS2 |
| `ds1norm.png` | Normalized data DS1 |
| `ds2norm.png` | Normalized data DS2 |

### **Batch Correction & Integration**

| File | Description |
|------|-------------|
| `afterbatcheffects.png` | UMAP after Harmony batch correction |
| `afterbatcheffectswithclusters.png` | Post-integration with clusters labeled |
| `nobatcheffects.png` | UMAP before batch correction (showing batch effect) |

### **Dimensionality Reduction**

| File | Description |
|------|-------------|
| `elbowplot.png` | PCA elbow plot (selecting PCs) |
| `PCheatmap1.png` | Heatmap of top PCs (genes driving variation) |
| `pcheatmap2.png` | Additional PC heatmap |

### **Clustering & Cell Type Identification**

| File | Description |
|------|-------------|
| `UMAPwithcluster.png` | UMAP with all 14 clusters labeled |
| `umapaftercleaning.png` | UMAP after contamination removal |
| `seuratcluster.png` | Seurat clustering results |
| `seuratumap.png` | UMAP colored by Seurat clusters |
| `tsnemap.png` | t-SNE visualization of clusters |
| `tsnemoreclusters.png` | t-SNE with additional resolution |
| `umapmoreclusters.png` | UMAP at higher clustering resolution |
| `umapwithnames.png` | UMAP with cell type names |

### **Marker Gene Expression**

| File | Description |
|------|-------------|
| `umapmarkergenes.png` | Key marker genes on UMAP (SOX2, EOMES, NEUROD2, etc.) |
| `tsnemarkergenes.png` | Marker genes on t-SNE |
| `withmarkergenesaft.png` | Markers after cleaning |
| `top10clustersheatmap.png` | Heatmap of top 10 markers per cluster |
| `heatmap3.png` | Additional marker heatmap |

### **Regional Identity (Dorsal vs Ventral)**

| File | Description |
|------|-------------|
| `dorsal vs ventral.png` | UMAP colored by region (Dorsal/Ventral) |
| `dvvumap.png` | Dorsal-ventral regional annotation |
| `degdvv.png` | Differential genes dorsal vs ventral |
| `deg2volcano.png` | Volcano plot: Dorsal vs Ventral DEGs |

### **Maturation Stages**

| File | Description |
|------|-------------|
| `maturationstage1uamp.png` | UMAP colored by maturation stage |
| `ipneuronprogenotorumap.png` | Progenitor/IP/Neuron stages on UMAP |
| `progenitor vs neuron vs cycling.png` | Three-way comparison |

### **Developmental Trajectory (Pseudotime)**

| File | Description |
|------|-------------|
| `psudotime.png` | Pseudotime values on UMAP |
| `spesudotime1.png` | Pseudotime visualization (subset) |
| `developmentalmarkeralongtrajectory.png` | Gene dynamics: GLI3, EOMES, NEUROD6 |
| `cellcycleremove.png` | Cell cycle regression effect |

### **Differential Gene Expression**

| File | Description |
|------|-------------|
| `degvolcanopvn.png` | Volcano: Progenitors vs Neurons |
| `topdegpvn.png` | Top DEGs: Progenitors vs Neurons |

### **Validation & Reference Comparison**

| File | Description |
|------|-------------|
| `celllevelslusters.png` | Cell-level validation clusters |
| `cellleveltranscriptome.png` | Transcriptome similarity (cell-level) |
| `umaptranscriptomelevel.png` | UMAP colored by reference correlation |

### **Gene Set Enrichment Analysis (GSEA)**

#### Progenitor Pathways
| File | Description |
|------|-------------|
| `progenotorenriched pathway1.png` | Top progenitor-enriched GO terms |
| `topprogenitorpathway.png` | Progenitor pathways (bar plot) |

#### Neuron Pathways
| File | Description |
|------|-------------|
| `neuronenrichedpathway.png` | Top neuron-enriched GO terms |
| `topenuronpathways.png` | Neuron pathways (bar plot) |

#### Regional Pathways
| File | Description |
|------|-------------|
| `dorsalpathway.png` | Dorsal lineage GO enrichment |
| `ventralpathway.png` | Ventral lineage GO enrichment |

#### Stage-Specific Pathways
| File | Description |
|------|-------------|
| `ipsecficipathways.png` | IP-specific pathways |
| `GOacrosspathways.png` | GO terms across all comparisons |
| `Geacrossstages.png` | Gene enrichment across maturation stages |
```

### Results - Tables (CSV)
```
results/tables/
├── deg_progenitors_vs_neurons.csv           # Progenitor vs Mature Neuron
├── deg_dorsal_vs_ventral.csv                # Dorsal vs Ventral lineages
├── deg_progenitor_vs_ip.csv                 # Progenitor vs Intermediate Progenitor
├── deg_ip_vs_neuron.csv                     # Intermediate Progenitor vs Neuron
├── deg_progenitor_vs_neuron_trajectory.csv  # Full trajectory comparison
├── gsea_progenitor_pathways.csv             # Progenitor GO enrichment
├── gsea_neuron_pathways.csv                 # Neuron GO enrichment
├── gsea_dorsal_pathways.csv                 # Dorsal lineage pathways
├── gsea_ventral_pathways.csv                # Ventral lineage pathways
├── gsea_ip_specific_pathways.csv            # Intermediate Progenitor pathways
└── pseudotime_results.csv                   # Pseudotime values and gene dynamics
```



### Python Export
```
exports/
└── integrated_clean.h5ad          # AnnData format for scVelo/scanpy
```

## How to Use This Analysis

### Running the Analysis

Execute the R markdown file:

```R
# Install packages first (one-time)
# Then render entire document:
rmarkdown::render("scrna_brain.Rmd", output_format = "html_document")

# Or run chunk-by-chunk interactively:
# - Open scrna_brain.Rmd in RStudio
# - Click "Run" on each code chunk
# - Explore outputs between runs
```



``
```
