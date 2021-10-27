# -*- coding: utf-8 -*-
"""Dimensional reduction and batch correction using Harmony."""
if __name__ == "__main__":
    import anndata as ad
    import matplotlib.pyplot as plt
    import scanpy as sc
    import scanpy.external as sce
    from helpers.logs.get_logger import get_logger
    from helpers.logs.sc_logs import set_sc_log
    from helpers.select_pcs import select_pcs

    LOG = snakemake.log[0]  # noqa: F821
    PARAMS = snakemake.params  # noqa: F821
    INPUT = snakemake.input  # noqa: F821
    OUTPUT = snakemake.output  # noqa: F821
    THREADS = snakemake.threads  # noqa: F821

    logger = get_logger(__name__, LOG)
    sc.settings = set_sc_log(sc.settings, logfile=LOG)
    sc.settings.n_jobs = THREADS

    # Concatenate samples
    adata = ad.concat(
        [sc.read_h5ad(path) for path in INPUT["data"]],
        join="outer",
        merge="same",
        label=None,
    )
    adata.obs_names_make_unique()
    logger.info(f"Adata read from {INPUT['data']}")
    logger.info(f"Input data: {adata}")

    # HVGs
    # Before normalisation as seurat_v3 expects raw counts
    if not PARAMS["nHVG"]:
        nHVG = max(min(len(adata.obs) / 2, 10000), 1000)
        logger.info(f"nHVG not provided. Using {nHVG}.")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=nHVG,
        flavor="seurat_v3",
        batch_key="lane",
        subset=False,
    )
    _ = sc.pl.highly_variable_genes(
        adata,
        log=False,
        show=False,
        save=False,
    )
    plt.savefig(OUTPUT["hvg"], dpi=300, bbox_inches="tight")
    plt.close()

    # Normalise
    # Exclude highly expressed to prevent skew of normalisation
    sc.pp.normalize_total(adata, exclude_highly_expressed=True)
    sc.pp.log1p(adata)

    # Save raw and filter
    adata.raw = adata

    # Regress and scale
    # No batch - covered with bbknn
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], n_jobs=None)
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    _ = sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False, save=False)
    plt.savefig(OUTPUT["elbow"], dpi=300, bbox_inches="tight")
    plt.close()

    # Harmony for batch correction
    # As it runs on all pcs include, we must first filter to desired
    npc, adata = select_pcs(adata, threshold=PARAMS["var_thresh"])
    logger.info(f"{npc} PCs used.")

    sce.pp.harmony_integrate(
        adata,
        key="lane",
        adjusted_basis="X_harmony",
        max_iter_harmony=50,
    )

    # And save
    adata.write_h5ad(OUTPUT["data"])
