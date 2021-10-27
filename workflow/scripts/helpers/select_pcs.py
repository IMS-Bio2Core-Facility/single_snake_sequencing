# -*- coding: utf-8 -*-
"""A semi-quantitative method for selecting nPC."""
from typing import Tuple

import numpy as np
from anndata._core.anndata import AnnData


def select_pcs(adata: AnnData, threshold: float = 0.75) -> Tuple[int, AnnData]:
    """Select the number of PCs using a cumulative variance threshold.

    Selecting the number of PCs to be kept in downstream scRNAseq analyses is not yet
        a solved problem. Fortunately, evidence suggests the analysis is _reasonably_
        robust to the number selected, so a semi-quantitative manner suffices. Here,
        I find the cumulative variance explained by the PCs and keep those PCs that
        explain up to a threshold of the variance. I find 0.75 to be a reasonable
        threshold if 40 PCs are initially calculated, though inspection of a knee
        plot can quide this decision.

    Note
    ----
    Currently, the function assumes that the PCA has been stored to the standard slots
        within the AnnData object. That is, the loadings are at `adata.varm["PCs"]`,
        the coordinates are at `adata.obsm["X_pca"]`, and the metadata at
        `adata.uns["pca"]`. If you call `sc.pl.pca` with the default parameters,
        this **will** be the case.

    Parameters
    ----------
    adata: AnnData
        An AnnData object containing a PCA calculation.
    threshold: float
        The cumulative variance threshold below which all PCs will be kept.

    Returns
    -------
    Tuple[int, AnnData]
        A tuple containing the number of selected PCs and the filtered AnnData object

    """
    var = adata.uns["pca"]["variance"]
    var_per = var / var.sum()
    npc: int = np.searchsorted(np.cumsum(var_per), threshold) + 1
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :npc]
    adata.varm["PCs"] = adata.varm["PCs"][:, :npc]
    return npc, adata
