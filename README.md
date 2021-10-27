# single_snake_sequencing - sc/snRNAseq Snakemake Workflow

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CI/CD](https://github.com/IMS-Bio2Core-Facility/single_snake_sequencing/actions/workflows/cicd.yml/badge.svg)](https://github.com/IMS-Bio2Core-Facility/single_snake_sequencing/actions/workflows/cicd.yml)
[![Codestyle: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Codestyle: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

A [Snakemake][sm] workflow for standardised sc/snRNAseq analysis.

Every single cell analysis is slightly different.
This represents what I would call a "core" analysis,
as nearly every analysis I perform start with something very akin to this.
Given this custom nature of single cell,
this workflow is not designed to be all encompassing.
Rather,
it aims to be **extensible**, **modular**, and **reproducible**.
Any given step can be easily modified -
as they are all self contained scripts -
and a new rule can be easily added -
see the [downstream rules](workflow/rules/downstream.smk) for an example.
Finally,
by taking advantage of the integrated [Conda][conda] and [Singularity][sing] support,
we can run the whole thing in an isolated environment.

## Notes on Installation

HOLDING

## Notes on Configuration

> :warning: Be sure to change the configuration to suit your project!

For a full discussion of configuration,
please see the [configuration README](config/README.md).

Briefly,
the general configuration file must be located at `config/config.yaml`.
A samplesheet containing information pertaining to the data must be supplied as well.
Both are schema validated.

## Notes on Data

This pipeline expects de-multiplexed fastq.gz files,
normally produced by some deriviative of `bcl2fastq` after sequencing.
They can (technically) be placed anywhere,
but we recommend creating a `data` directory in your project for them.

## Notes on the tools

The analysis pipeline was run using Snakemake v6.6.1.
The full version and software lists can be found under the relevant yaml files in `workflow/envs`.
The all reasonable efforts have been made to ensure that the repository adheres to the best practices
outlined [here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).

## Notes on the analysis

For a full discussion on the analysis methods,
please see the [technical documentation](workflow/documentation.md).

Briefly,
the count matrix was produced using Cellranger,
droplet calling with `DropletUtils::emptyDrops`,
doublet detection with `SOLO` from the `scVI` family,
batch effect removal with `harmonypy`,
and general analysis and data handling with `scanpy`.

## Future work

- [ ] Supply tests
- [ ] Track lane in samples that have been pooled and de-multiplexed
- [ ] Parallelise emptyDrops
- [ ] Support custom references
- [ ] Support SCTransform?

[sm]: https://snakemake.readthedocs.io/en/stable/index.html "Snakemake"
[conda]: https://docs.conda.io/en/latest/ "Conda"
[sing]: https://sylabs.io/singularity/ "Singularity"
