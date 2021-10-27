# Technical Documentation

What follows is a technical documentation of each step in the pipeline.

The steps `get_cellranger`, `get_reference`, `unpack_bcl2fastq`, and `clean_names` are not here covered,
as the simply move or rename files.

## Rule: counts

Uses [Cellranger counts][cr] to generate the counts matrix for each sample.

The `--noseconday` flag is passed,
since we are (obviously) performing our own.

The `--include-introns` flag is controlled by a user-provide boolean switch,
allowing flexibility for choosing between single-cell and single-nucleus data.

Currently,
only a `--sample` flag is passed,
collating all counts into a single matrix.
Future work will also pass a `--lane` flag,
which will allow for batch tracking,
even when samples have been pooled and de-multiplexed.

The `--expect-cells` flag takes an integer count that must be specified by the user
based on knowledge about experimental design.
In my experience,
Cellranger always over estimates the number of cells,
regardless of this flag,
hence the custom implementation of `emptyDrops`.

## Rule: filter_empty

Runs [DropletUtils::emptyDrops][empty] on the 10x output.

The challenge in implementing this function is determining the correct value for `lower`,
Set it too high,
and you will lose droplets that are actual cells and confound the identification of the rest.
Set it too low,
and you will gain cells at the cost of increased error in the cell call
(as garbage will be included).
In my experience,
the most effective way to set this in an unbiased manner is to take the median of counts
between 25 and 2000 (to prevent 0 inflation)
and set the bound at 2 standard deviations above this mean.

## Rule: qc

Standard sc/snRNAseq quality control metrics run using [Scanpy][scanpy].

This is run on a per sample basis.
Filtering occurs on the following criteria:

1. Any droplet whose FDR from `filter_empty` is > 0.005 is dropped.
1. Any feature not found in at least 1 cell is dropped.
1. Any cell with more than the user-defined threshold of `percent mitochondrial counts` is dropped.
  - This defaults to 5%. This is ~ correct for single nucleus data. 20% is better for single cell.
1. Any cell with more than the user-defined threshold of `total counts` is dropped
1. Any cell with less than the user-defined threshold of `genes by counts` is dropped
1. If Cellranger included introns, then all mitochondrial genes are dropped.

Finally,
the sample label is extracted from the Snakemake wildcard
to label the sample for future downstream analysis.

## Rule: doublets

Doublet detection performed using [Solo][solo] implemented in the [scVI][scvi] package.

Briefly,
Solo uses an (unsupervised) variational auto-encoder to embed cells
then appends a feed-forward neural network to form a supervised classifier.
Please see their original [manuscript][solo_paper] for more details.

As a neural network,
it has the usual litany of strengths and weaknesses.
It is GPU accelerable,
allowing it perform well on large datasets.
However,
it can be prone to over-fitting,
a problem addressed in this pipeline by the use of early-stopping.

Solo is run on a lane-by-lane basis,
as each lane will have different technical variation.

## Rule: dim_reduc

Dimensional reduction and bath-correction using [Scanpy][scanpy] and [Harmony][harmony].

At this point in the analysis,
all samples and lanes are concatenated into a single [Anndata][adata] object.
Lane and sample labels are retained for identification downstream,
where necessary.

The number of highly variable genes to use is a key choice in sc/snRNAseq analysis,
and one that frequently can feel quite arbitrary.
In effort to make this more quantitative,
I default to half the number of cells or 10,000 gene,
whichever is less,
while enforcing a floor of 1000 genes
Should the user desire,
this default can be overridden by specifying `nHVG` in the
[configuration file](config/config.yaml).

Highly expressed genes -
those accounting for more than 5% of counts in cells -
are excluded from the calculation of the normalisation factor
to prevent skew in the normalised values of other genes.
Frequently,
this simply excludes MALAT1,
though the user may check the log files for the exact details.

After a `log1p` transformation,
the effects of `total counts` and `percent counts from mitochondria` are regressed out.
The latter can be debated as to whether it is a relevant step.
I find that,
unless there is a strong reason to suspect a phenotype defined by mitochondrial content,
it reduces noise in the final data.
Additionally,
for single nucleus data,
these genes can be considered a source of technical error,
as they should not be detected.

50 PCs are calculated on the variable features.
To automate the number of PCs to keep,
those PCs accounting for 85% of the variance are kept.
This can be overridden by the user,
should there be a motivation to keep more or less PCs.

Batch correction is performed over lanes using [Harmony][harmony].
The correct application of batch correction is a discussion beyond the scope of these technical notes.
A number of papers (see [here][harmony_paper_1] and [here][harmony_paper_2])
suggest that Harmony is a robust choice for batch correction.
Here,
it is run on all PCs selected in the previous step until model convergence.

## Rule: cluster

Clustering using [Leiden][leiden] with [Scanpy][scanpy].

A neighbor network is calculated using 15 neighbors on the Harmony-corrected PCs.
The [Leiden algorithm][leiden] for clustering is used on this network
at a user-defined resolution.
A basic differential expression test is performed in a 1-vs-all manner across all clusters,
using a T-test with overestimated variance and a Benjamini-Hochberg multiple correction.
The user is advised that these should be consider with a pinch of salt,
as they are **not** batch corrected.
For more robust differential expression analysis across a condition or cell type,
the user is directed to tools such as [MAST][mast] that can be integrated into this pipeline.

## Rule: densities

Density embeddings using [Scanpy][scanpy].

This step allows the user to definge categorical covariates on which
density embeddings will be calculated.
I find this a useful sanity check on my analysis,
as well as providing an example of how the basic workflow may be extended.

[cr]: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count "Cellranger"
[empty]: https://rdrr.io/bioc/DropletUtils/man/emptyDrops.html "emptyDrops"
[scanpy]: https://scanpy.readthedocs.io/en/stable/ "Scanpy"
[solo]: https://github.com/calico/solo "Solo"
[scvi]: https://scvi-tools.org/ "scVI"
[solo_paper]: https://www.sciencedirect.com/science/article/pii/S2405471220301952 "Solo Manuscript"
[harmony]: https://github.com/slowkow/harmonypy "Harmony"
[adata]: https://anndata.readthedocs.io/en/latest/index.html "Anndata"
[harmony_paper_1]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9 "Tran, et al."
[harmony_paper_2]: https://academic.oup.com/nar/article/49/7/e42/6125660 "Chazarra-Gil, et al."
[leiden]: https://www.nature.com/articles/s41598-019-41695-z "Leiden clustering"
[mast]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5 "MAST"
