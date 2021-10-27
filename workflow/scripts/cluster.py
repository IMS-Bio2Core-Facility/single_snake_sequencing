# -*- coding: utf-8 -*-
"""Cluster single cell data."""
if __name__ == "__main__":
    from pathlib import Path

    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc
    from helpers.logs.get_logger import get_logger
    from helpers.logs.sc_logs import set_sc_log
    from helpers.write_de import get_de
    from matplotlib.backends.backend_pdf import PdfPages

    LOG = snakemake.log[0]  # noqa: F821
    PARAMS = snakemake.params  # noqa: F821
    INPUT = snakemake.input  # noqa: F821
    OUTPUT = snakemake.output  # noqa: F821
    THREADS = snakemake.threads  # noqa: F821

    logger = get_logger(__name__, LOG)
    sc.settings = set_sc_log(sc.settings, logfile=LOG)
    sc.settings.n_jobs = THREADS

    adata = sc.read_h5ad(INPUT["data"])
    logger.info(f"Adata read from {INPUT['data']}")
    logger.info(f"Input data: {adata}")

    sc.pp.neighbors(adata, use_rep="X_harmony")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=PARAMS["res"])
    logger.info(
        f"{len(adata.obs.leiden.unique())} clusters identified at resolution {PARAMS['res']}"
    )

    _ = sc.pl.umap(
        adata,
        color=[
            "lane",
            "sample",
            "leiden",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo_prot",
        ],
        use_raw=True,
        show=False,
        save=False,
        ncols=2,
    )
    plt.savefig(OUTPUT["umap"], dpi=300, bbox_inches="tight")
    plt.close()

    # Marker genes
    if len(PARAMS["markers"]) == 0:
        Path(OUTPUT["markers"]).touch(exist_ok=True)
        logger.info("No markers provided. File touched, but empty.")
    else:
        with PdfPages(OUTPUT["markers"]) as pdf:
            for i in range(0, len(PARAMS["markers"]), 9):
                genes = PARAMS["markers"][i : i + 9]
                genes = [x for x in genes if x in adata.var_names]
                fig = sc.pl.umap(
                    adata,
                    color=genes,
                    use_raw=True,
                    return_fig=True,
                    save=False,
                    ncols=3,
                )
                pdf.savefig(fig, dpi=300, bbox_inches="tight")
                plt.close()
                logger.info(f"Markers {genes} plotted.")

    # Run dendrogram
    sc.tl.dendrogram(
        adata,
        use_rep="X_harmony",
        groupby="leiden",
        optimal_ordering=True,
    )

    # Differential expression
    # Not ranking by abs will return the 100 highest scoring upregulated genes
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        use_raw=True,
        method="t-test_overestim_var",
        pts=True,
    )

    # And plot
    _ = sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=5,
        values_to_plot="logfoldchanges",
        save=False,
        show=False,
        cmap="coolwarm",
        vmin=-4,
        vmax=4,
        colorbar_title="log2FC",
    )
    plt.savefig(OUTPUT["dot"], dpi=300, bbox_inches="tight")
    plt.close()

    # Write results
    results = get_de(adata)
    with pd.ExcelWriter(OUTPUT["de"], engine="openpyxl") as writer:
        for j, df in results:
            df.to_excel(writer, sheet_name=f"cluster_{j}")

    logger.info(f"Final adata: {adata}")
    adata.write_h5ad(OUTPUT["data"])
