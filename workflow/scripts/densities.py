# -*- coding: utf-8 -*-
"""Calculate embedding densities for sc/snRNAseq data."""
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import scanpy as sc
    from helpers.logs.get_logger import get_logger
    from helpers.logs.sc_logs import set_sc_log
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

    # Calculate
    with PdfPages(OUTPUT["densities"]) as pdf:
        for feature in PARAMS["features"]:
            sc.tl.embedding_density(adata, groupby=feature)
            fig = sc.pl.embedding_density(
                adata,
                groupby=feature,
                group="all",
                color_map="plasma",
                ncols=3,
                save=False,
                return_fig=True,
            )
            fig.suptitle(f"Densities for feature: {feature}")
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close()
            logger.info(f"Feature {feature} plotted")

    # Save data
    logger.info(f"Final data: {adata}")
    adata.write_h5ad(OUTPUT["data"])
