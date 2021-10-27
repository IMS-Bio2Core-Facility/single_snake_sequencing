# -*- coding: utf-8 -*-
"""Nox session configuration."""
import nox
from nox.sessions import Session

PACKAGE: str = "workflow/scripts"
LOCATIONS: list[str] = [
    PACKAGE,
    "noxfile.py",
]
VERSIONS: list[str] = ["3.9"]
PIP_PARAMS: list[str] = ["-r"]

nox.options.stop_on_first_error = True
nox.options.default_venv_backend = "virtualenv"
nox.options.reuse_existing_virtualenvs = True
nox.options.sessions = [
    "form",
    "lint",
    "type",
    "security",
]


@nox.session(python=VERSIONS[0])
def form(session: Session) -> None:
    """Format code with isort and black."""
    args = session.posargs or LOCATIONS
    session.install(*PIP_PARAMS, "workflow/envs/form.txt")
    session.run("isort", *args)
    session.run("black", *args)


@nox.session(python=VERSIONS)
def lint(session: Session) -> None:
    """Lint files with flake8."""
    args = session.posargs or LOCATIONS
    session.install(*PIP_PARAMS, "workflow/envs/lint.txt")
    session.run("pflake8", *args)


@nox.session(python=VERSIONS)
def type(session: Session) -> None:
    """Type check files with mypy."""
    args = session.posargs or LOCATIONS
    session.install(*PIP_PARAMS, "workflow/envs/type.txt")
    session.run(
        "mypy",
        "--ignore-missing-imports",
        "--disable-error-code",
        "name-defined",
        *args
    )


# The snakemake install is complex, and runs better in conda
@nox.session(python=VERSIONS)
def security(session: Session) -> None:
    """Check security safety."""
    args = session.posargs or []
    session.install("safety")
    session.run("safety", "check", "-r", "workflow/envs/security.txt", *args)
