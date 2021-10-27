# -*- coding: utf-8 -*-
"""Set logging properties for scanpy."""
from typing import Literal

import scanpy as sc


def set_sc_log(
    settings: sc._settings.ScanpyConfig,
    logfile: str = "log.txt",
    verbosity: Literal[0, 1, 2, 3, 4] = 4,
) -> sc._settings.ScanpyConfig:
    """Set logging properties for scanpy.

    Parameters
    ----------
    settings : scanpy._settings.ScanpyConfig
        The scanpy settings instance to configure.
    logfile : str
        Where to log results to
    verbosity : Literal[0,1,2,3,4]
        How verbose to be

    Returns
    -------
    sc._settings.ScanpyConfig
        The configured scanpy settings instance
    """
    settings.logfile = logfile
    settings.verbosity = verbosity
    return settings
