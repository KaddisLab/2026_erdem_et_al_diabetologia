# Input data

The pipeline needs one file in this directory:

| | |
|---|---|
| Filename | `060425_scRNA_v3.3.rds` |
| Size | ~9 GB |
| Format | RDS containing a Seurat v5 object |
| Contents | PanKBase single-cell RNA-seq, ~448,000 cells from 139 donors |

## Download

The PanKBase Seurat object for the publicly available scRNA-seq datasets is
archived at <https://pubrepo.coh.org/publication/17>.

Place the file here, so that `data/060425_scRNA_v3.3.rds` exists relative to
the repository root, then run `Rscript reproduce_figures.R` from that root.

## Check the file loaded correctly

```r
obj <- readRDS("data/060425_scRNA_v3.3.rds")
ncol(obj)                     # ~448,000 cells
length(unique(obj$aliases))   # 139 donors
```

Loading the object needs roughly 20 GB of RAM; the UCell scoring step that
follows is the memory peak for the pipeline as a whole.

## Provenance

The object integrates data from three sources. Only the first two contribute to
the analysis cohort:

- **HPAP** — Human Pancreas Analysis Program, <https://hpap.pmacs.upenn.edu/>
- **IIDP** — Integrated Islet Distribution Program, <https://iidp.coh.org/>
- **Prodo Laboratories** — cytokine-treated cells, excluded here by the
  control/vehicle treatment filter

`filter_pankbase_main()` keeps only control and vehicle treatments with a
curated diabetes status, and `filter_pankbase_ductal()` then keeps ductal cell
types, giving the cohort reported in the paper: **20,605 ductal cells from 86
donors**.

The bulk mRNA-seq data from the same study are deposited separately in GEO
under accession [GSE309294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309294).
