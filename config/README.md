# Configuration

The configuration keys that are expected are given below.
Don't worry about typos, etc.
These are all enforced with Snakemake's brilliant
[schema validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation).

## congif.yaml

### samplesheet

Path to the samplesheet.
This defaults to `config/samples.yaml`.

### get_cellranger

- url: str, required. Url to retrieve Cellranger

### get_reference

- url: str, required. Url to retrieve Cellranger reference

### counts

- introns: bool, required. Whether introns should be included (ie. sn- vs sc-RNAseq)
- n_cells: int, required. The expected number of cells to recover per sample
- mem: int, required. The memory, in Gb, available to **each instance** of cell ranger

### filter_empty

- niters: int, required. The number of iterations to run `DropletUtils::emptyDrops`

### qc

- pct_counts_mt: int, required. Cells with more than this value will be discarded
- total_counts: int, required. Cells with more than this value will be discarded
- n_genes_by_counts: int, required. Cell with less than this value will be discarded

### dim_reduc

This top level key is optional.

- nHVGs: int, optional. The number of variable genes to use. Defaults to 1/2 the number of cells or 10,000, whichever is less, but will never be below 1000.
- var_thresh: float, optional. Keep PCs explaining upto this fraction of variance. Between 0 and 1, defaults to 0.85.

### cluster

- res: float, required. The resolution to use for Leiden clustering.
- markers: list of strings, optional. Marker to plot expression for. No default.

### densities

This top level key is optional.

- features: list of strins, optional. Categorical features to calculate embeddings for. Defaults to "sample".

## samples.yaml

The top level keys are the lanes from the sequencer.
The second level keys are the samples from that lane.
The third level keys are the paths to the R1 and R2 data.
This schema also validated.
