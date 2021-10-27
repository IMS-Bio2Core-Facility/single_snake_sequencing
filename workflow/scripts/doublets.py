# -*- coding: utf-8 -*-
"""Run Solo with SCVI."""
if __name__ == "__main__":
    import sys

    import anndata as ad
    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc
    import scvi
    from helpers.logs.get_logger import get_logger
    from pandas.plotting import table

    LOG = snakemake.log[0]  # noqa: F821
    PARAMS = snakemake.params  # noqa: F821
    INPUT = snakemake.input  # noqa: F821
    OUTPUT = snakemake.output  # noqa: F821
    WILDS = snakemake.wildcards  # noqa: F821

    logger = get_logger(__name__, LOG)
    scvi.settings.seed = 42

    # Get data
    adata = ad.concat(
        [sc.read_h5ad(path) for path in INPUT["data"]],
        join="outer",
        merge="same",
        label=None,
    )
    adata.obs_names_make_unique()
    adata.obs["lane"] = str(WILDS["lane"])
    logger.info(f"Adata read from {INPUT['data']}")
    logger.info(f"Input data: {adata}")

    # SCVI doesn't seem to have a way to specify where the output logs to,
    # So we have to wrap and redirect
    with open(LOG, "a") as file:
        sys.stderr = sys.stdout = file

        # No covariates - not used with SOLO
        scvi.data.setup_anndata(
            adata,
        )
        vae = scvi.model.SCVI(adata)
        vae.train(
            early_stopping=True, check_val_every_n_epoch=2, early_stopping_patience=20
        )

        solo = scvi.external.SOLO.from_scvi_model(
            vae,
        )
        solo.train(train_size=0.9, early_stopping=True, early_stopping_patience=20)
        preds = solo.predict(soft=True)

    # SCVI appends a random -0 to the index...
    adata.obs.index = adata.obs.index + "-0"
    adata.obs["singlet_score"] = preds.loc[:, "singlet"]
    adata.obs["doublet_score"] = preds.loc[:, "doublet"]

    # Save a results table
    preds = (
        pd.DataFrame.from_dict(
            {
                "Sample": adata.obs["sample"],
                "Singlet": adata.obs["singlet_score"] > adata.obs["doublet_score"],
                "Doublet": adata.obs["singlet_score"] < adata.obs["doublet_score"],
            }
        )
        .groupby("Sample")
        .sum()
    )

    ax = plt.subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    _ = table(ax, preds, cellLoc="center")
    plt.savefig(OUTPUT["table"], dpi=300, bbox_inches="tight")
    plt.close()

    # Filter
    adata = adata[adata.obs.singlet_score > adata.obs.doublet_score, :]
    logger.info(f"Data after doublet calls: {adata}")

    # Clean up after scvi
    del adata.uns
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith("_")].copy()
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    adata.obs["lane"] = adata.obs["lane"].astype("category")

    adata.write_h5ad(OUTPUT["data"])
    logger.info(f"Data saved to {OUTPUT['data']}")
