## Multiomics atlas of mouse embryogenesis reveals coordinated germ layer development and critical transition at E10.5-E11.5
---
We generated a paired single-cell dual-omics atlas of mouse embryogenesis using HT-scCAT-seq, capturing joint transcriptome and chromatin-accessibility profiles for 100,999 cells across five stages (E8.5–E12.5) and annotating 33 cell types spanning all three germ layers and extraembryonic tissues.
![overview](https://github.com/vaporer/Multiomics-atlas-of-mouse-embryogenesis/blob/master/%E8%B5%84%E6%BA%90%201.png)


## Reproducing the results from the paper

```
00 Mouse embryo QC
01 Annotations verification
01 Multiomics integration
02 Stage related analysis
03 peak to gene linkage
04 Candidate enhencer
04 Creating custom cistarget database
04 Run pycistarget 
04 SCENICplus downstream analysis
05 Neuron analysis
```

## Softwares/Packages needed
**Python packages**
- numpy (v1.19.5)
- scipy (v1.5.3)
- pandas (v1.1.5)
- seaborn (v0.11.0.rc0)
- scanpy (v1.9.0)
- scikit-learn (v0.22.1)

**R packages**
- Matrix (v1.3-3)
- Signac (v1.6.0)
- Seurat (v4.1.0)
- monocle (v2.20.0)
- GenomeInfoDb (v1.30.1)
- EnsDb.Mmusculus.v79 (v2.99.0)
- EnsDb.Hsapiens.v86 (v2.99.0)
- ggplot2 (v3.35)
- patchwork (v1.1.1)
- dplyr (v1.0.8)
- stringr (v1.4.0.9000)
- stringi (v1.7.6)
- FigR (v0.3.5)
- scCustomize (v1.1.3)
- JASPAR2020 (v0.99.10)
- TFBSTools (v1.40.0)
- FNN (v1.1.4.1)
- cisTopic (v2.1.0)
- GenomicRanges (v1.54.1)
- SummarizedExperiment (v1.32.0)
- tidyverse (v2.0.0)
- igraph (2.0.3)
- ComplexHeatmap (v2.18.0)
- monocle3 (v1.3.1)
- cicero (v1.26.0)
- motifmatchr (v1.24.0)
- chromVAR (v1.24.0)

**other Bioinformatics packages**
- chromap (0.2.1-r369)
- MACS2 (v2.2.7.1)
- STARsolo (v2.7.9a)
- Cellranger-ATAC (v2.0.0)
- samtools (v1.9)
- bedtools (v2.30.0)
- seqtk (v1.3-r106)
- sratoolkit (v3.0.0)
- faSize, bedClip, calc, addCols from [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)
## Citation

 
=======
---
