# -*- coding: utf-8 -*-
"""Calculate QC metrics and filter."""
if __name__ == """__main__""":
    from pathlib import Path

    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc
    from helpers.logs.get_logger import get_logger
    from helpers.logs.sc_logs import set_sc_log

    LOG = snakemake.log[0]  # noqa: F821
    PARAMS = snakemake.params  # noqa: F821
    INPUT = snakemake.input  # noqa: F821
    OUTPUT = snakemake.output  # noqa: F821
    WILDS = snakemake.wildcards  # noqa: F821

    logger = get_logger(__name__, LOG)
    sc.settings = set_sc_log(sc.settings, logfile=LOG)

    adata = sc.read_10x_mtx(INPUT["dir"], make_unique=True)
    adata.obs_names_make_unique()
    logger.info(f"Adata read from {INPUT['dir']}")
    logger.info(f"Input data: {adata}")

    empty = pd.read_csv(INPUT["empty"], index_col=0)
    empty.index = adata.obs_names
    empty.FDR = empty.FDR.fillna(1)
    logger.info(f"Empty read from {INPUT['empty']}")

    # Transfer call
    adata.obs["is_cell"] = empty.FDR < 0.005

    # Initial filtering
    adata = adata[adata.obs["is_cell"], :]
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)

    # Mark mitochondrial and ribosomal genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo_prot"] = adata.var_names.str.startswith("RPS", "RPL")

    # Calculate qc metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo_prot"], percent_top=None, log1p=False, inplace=True
    )

    # Various plots
    # Specify various sc figure settings
    sc.set_figure_params(dpi_save=300, fontsize=12, figsize=(6, 6), format="png")
    sc.settings.figdir = Path(".")

    _ = sc.pl.violin(
        adata,
        [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ],
        multi_panel=True,
        show=False,
        save=False,
    )
    plt.savefig(OUTPUT["qc_violin"], dpi=300, bbox_inches="tight")
    plt.close()

    _ = sc.pl.scatter(
        adata,
        "total_counts",
        "pct_counts_mt",
        show=False,
        save=False,
    )
    plt.savefig(OUTPUT["mt_scatter"], dpi=300, bbox_inches="tight")
    plt.close()

    _ = sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        show=False,
        save=False,
    )
    plt.savefig(OUTPUT["n_genes"], dpi=300, bbox_inches="tight")
    plt.close()

    # Filter cells
    adata = adata[adata.obs.pct_counts_mt < PARAMS["pct_counts_mt"], :]
    adata = adata[adata.obs.total_counts <= PARAMS["total_counts"], :]
    adata = adata[
        adata.obs.n_genes_by_counts >= PARAMS["n_genes_by_counts"],
        :,
    ]

    # Filter genes
    # If introns, then single nucleus -> mt is artefact, so drop
    if PARAMS["introns"]:
        adata = adata[:, ~adata.var.mt]

    # Final plot...
    _ = sc.pl.highest_expr_genes(adata, n_top=50, show=False, save=False)
    plt.savefig(OUTPUT["n_top"], dpi=300, bbox_inches="tight")
    plt.close()

    # Mark and save
    adata.obs["sample"] = WILDS["sample"]
    adata.write_h5ad(OUTPUT["data"])
    logger.info(f"Filtered data: {adata}")
    logger.info(f"Data saved to {OUTPUT['data']}")
