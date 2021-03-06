# single_snake_sequencing - sc/snRNAseq Snakemake Workflow

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CI/CD](https://github.com/IMS-Bio2Core-Facility/single_snake_sequencing/actions/workflows/cicd.yaml/badge.svg)](https://github.com/IMS-Bio2Core-Facility/single_snake_sequencing/actions/workflows/cicd.yaml)
[![Codestyle: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Codestyle: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

A [Snakemake][sm] workflow for standardised sc/snRNAseq analysis.

**_Find us in the "Standardized Usage" Section of the [Snakemake Workflow Catalog][sm_wc]_**

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

A full walkthrough on how to install and use this pipeline
can be found [here][sm_wc].

To take advantage of Singularity,
you'll need to install that separately.
If you are running on a Linux system,
then singularity can be installed from conda like so:

```shell
conda install -n snakemake -c conda-forge singularity
```

It's a bit more challenging for other operating systems.
Your best bet is to follow their instructions
[here][sing_install].
But don't worry!
**Singularity is _not_ regquired!**
Snakemake will still run each step in its own Conda environment,
it just won't put each Conda environment in a container.

### Get the Source Code

Alternatively,
you may grab the source code.
You likely will not need these steps if you aren't planning to contribute.

Navigate to our [release][releases]
page on github and download the most recent version.
The following will do the trick:

```shell
curl -s https://api.github.com/repos/IMS-Bio2Core-Facility/single_snake_sequencing/releases/latest |
grep tarball_url |
cut -d " " -f 4 |
tr -d '",' |
xargs -n1 curl -sL |
tar xzf -
```

After querying the github api to get the most recent release information,
we grep for the desired URL,
split the line and extract the field,
trim superfluous characters,
use `xargs` to pipe this to `curl` while allowing for re-directs,
and un-tar the files.
Easy!

Alternatively,
for the bleeding edge,
please clone the repo like so:

```shell
git clone https://github.com/IMS-Bio2Core-Facility/single_snake_sequencing
```

> :warning: **Heads Up!**
> The bleeding edge may not be stable,
> as it contains all active development.

However you choose to install it,
`cd` into the directory.

## Notes on Data

This pipeline expects de-multiplexed fastq.gz files,
normally produced by some deriviative of `bcl2fastq` after sequencing.
They can (technically) be placed anywhere,
but we recommend creating a `data` directory in your project for them.

## Notes on the tools

The analysis pipeline was run using Snakemake v6.6.1.
The full version and software lists can be found under the relevant yaml files in `workflow/envs`.
The all reasonable efforts have been made to ensure that the repository adheres to the best practices
outlined [here][sm_bp].

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
[mambaforge]: https://github.com/conda-forge/miniforge#mambaforge "Mambaforge"
[sing_install]: https://sylabs.io/guides/3.8/admin-guide/installation.html#installation-on-windows-or-mac "Singularity Install"
[sm_wc]: https://snakemake.github.io/snakemake-workflow-catalog?usage=IMS-Bio2Core-Facility/single_snake_sequencing "Usage Instructions"
[sm_bp]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html "Snakemake Best Practices"
[releases]: https://github.com/IMS-Bio2Core-Facility/single_snake_sequencing/releases "Releases"
