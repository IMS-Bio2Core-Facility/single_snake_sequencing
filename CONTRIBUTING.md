# Contributing

> :warning: **TLDR**
>
> Create a conda environment with:
>
>     ```shell
>     conda env create -f workflows/envs/devel.yml
>     conda activate devel
>     pre-commit install --install-hooks
>     ```
>
> then commit with `cz commit`,
> and everything should take care of itself!

> :warning: Linux-Dependent
>
> Due to dependencies on Singularity and Nvidia,
> these steps are Linux-dependent.
> A workflows/envs/devel_macos.yaml is provided,
> but it will **not** support Nvidia/Cuda/GPU acceleration.

Welcome, friend!
Open-source software isn't open open-source without the community.
We appreciate your interest and welcome all contributions.
While your here,
we respectfully ask that you abide by our [code of conduct](./CODE_OF_CONDUCT.md).
To help keep everything moving smoothly,
we have a few guidelines.

## Bugs

If you think you've found a bug,
let us know [here][issues].
We'll do our best to deal with it ASAP,
but please be patient as we also work many other projects!

## Developing

If you think you can fix one of the bugs,
or would like to submit a new feature,
then let's get coding!

Once you've cloned the repository,
fork it,
and get your development environment set up.
We use [conda][conda],
[nox][nox],
and [pre-commit][pre-commit]
to handle environments, testing, and linting.
Tests are handled by [snakemake][snakemake].
Between them,
they make sure that all checks run in isolated environments.
Please make sure you activate them before making any changes,
using the commands below:

```shell
conda env create -f workflows/envs/devel.yml
conda activate bic083-devel
pre-commit install --install-hooks
```

Now,
you don't even have to think about linting or testing.
When you commit a change,
pre-commit will automatically run [black][black],
[isort][isort],
[mypy][mypy],
a suite of [flake8][flake8]-based linters,
[snakefmt][snakefmt],
and `snakemake --lint`.
When you push a change,
github actions will trigger a more robust suit using Nox,
including security check and snakemake testing.

Sometimes,
its useful to run these lints manually.
The easiest way to do this is to use pre-commit:

```shell
pre-commit run -a
```

Which will run all lints except the slower security and test checks.
To run these,
use:

```shell
nox -s security
```

for the former and:

```
pytest .tests/unit # for unit tests
pytest .tests/integration # for integration tests
```

for the latter.

### Commits

We use a GitHub action for
[python-semantic-release][psr]
to manage our version numbers.
This automatically parses your commit messages to determine if a new release is nesessary,
and, if so, what kind (ie. major, minor, patch).
So your commit message is very important!

But we also don't want you stressing about how to format your commit message.
To that end,
we use a python implementation of
[commitizen][cz]
to handle that for you!
Just commit using `cz commit` instead of `git commit`,
and enjoy the magic!

### Documentation and Testing

Speaking of documentation and testing -
if you add new code,
please add documentation and tests for it as well.
We use [napoleon numpy][docstrings]
for our docstrings.

Testing a workflow is a bit different than testing a software package,
though the same principles apply.
To generate unit tests after changing the analysis scripts,
first run the workflow on a small test dataset.
Then run:

```shell
snakemake --generate-unit-tests
```

Be careful,
as this will reset all existing unit tests.
These can be recovered with `git reset`,
so don't worry.
Then,
using the existing tests as a guide,
tweak your new tests until they run.

For integration testing,
include all the shell code,
data,
configuration files,
etc.,
necessary for a minimal run under `.tests/integration`.

For further details,
see Snakemake's own guides [here][unit_tests] and [here][repro].

## Review

Once your happy with your code,
open a pull-request,
and we will reveiw ASAP.
If you pull-request is not passing on github actions,
or something else confuses us
(we are, after all, only human!),
we might ask for some small changes.
Once everything is looking good,
we will merges in your changes!

## From the Command-Line

> :fireworks: **Optional fun!**
> None of this section is required, but we find it useful and hope you do, too!

If you haven't heard of it already,
give a peek to to [gh][gh],
GitHub's official CLI.
It allows to manage all of the above steps from the command-line,
from forking,
to raising issues,
and checking on the status of your pull request.
Not a necessity,
but for you terminal warriors out there,
it just might help!

[issues]: https://github.com/IMS-Bio2Core-Facility/BIC089-Stefania-Carrobio/issues "Issues"
[conda]: https://docs.conda.io/en/latest/ "Conda"
[nox]: https://nox.thea.codes/en/stable/ "Nox"
[pre-commit]: https://pre-commit.com/ "Pre-commit"
[snakemake]: https://snakemake.readthedocs.io/en/stable/index.html "Snakemake"
[black]: https://github.com/psf/black "Black"
[isort]: https://pycqa.github.io/isort/ "iSort"
[mypy]: https://mypy.readthedocs.io/en/stable/index.html "Mypy"
[flake8]: https://flake8.pycqa.org/en/latest/ "Flake8"
[snakefmt]: https://github.com/snakemake/snakefmt#github-actions "snakefmt"
[psr]: https://github.com/relekang/python-semantic-release "Python Semantic Release"
[cz]: https://commitizen-tools.github.io/commitizen/index.html "Commitizen"
[docstrings]: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html "Numpy Napoleon Docstrings"
[unit_tests]: https://snakemake.readthedocs.io/en/stable/snakefiles/testing.html#snakefiles-testing "Snakemake Unit Tests"
[repro]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html "Snakemake Reproducibility"
[gh]: https://github.com/cli/cli "GitHub CLI"
