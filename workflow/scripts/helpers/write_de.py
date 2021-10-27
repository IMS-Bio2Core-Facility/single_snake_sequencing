# -*- coding: utf-8 -*-
"""Functions for processing differential expression data."""
from typing import List, Tuple

import pandas as pd
from anndata._core.anndata import AnnData


def get_de(
    adata: AnnData, key: str = "rank_genes_groups", group_by: str = "leiden"
) -> List[Tuple[str, pd.core.frame.DataFrame]]:
    """Get DE data from adata.uns as pd.Dataframes.

    While it stores alot of information, the adata.uns['rank_gene_groups']
    slot is un-nvaigable and un-readable in its current form. This is a helper
    function to get into dataframes for easy manipulation or exporting

    Parameters
    ----------
    adata: AnnData
        The AnnData object containing the DE test results
    key: str
        Default: "rank_gene_groups"
        Where the results are stored in adata.uns
    group_by: str
        Default: "Leiden"
        The clustering used to generate the data. This used to retrieve
        the cluster names

    Returns
    -------
    List[Tuple[pd.core.frame.DataFrame]]
        A list of tuples containing the cluster name and results.
    """
    # Get requisite initial data
    clusters = sorted(adata.obs[group_by].unique())
    de = adata.uns[key]
    enrichment = de["pts"] / de["pts_rest"]

    # Get DE data
    dfs = [
        (
            x,
            pd.DataFrame(
                {
                    "Gene": de["names"][x],
                    "Z-score": de["scores"][x],
                    "Adj. p_val": de["pvals_adj"][x],
                    "Log2 FC": de["logfoldchanges"][x],
                }
            )
            .set_index("Gene")
            .pipe(lambda df: df[df["Z-score"].ge(0)])
            .join(enrichment.loc[:, x])
            .rename(columns={x: "Enrichment"})
            .head(100),
        )
        for x in clusters
    ]

    return dfs
