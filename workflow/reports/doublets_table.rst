Counts the number of singlets and doublets per sample. Doublets were determined using the scVI_ implementation of SOLO_. Briefly, SOLO use a variational autoencoder to embed cells *unsupervised* before appending a neural network to form a *supervised* classifier. This classifier is then trained on simulated doublet created from the dataset.

.. _scVI: https://scvi-tools.org/
.. _SOLO: https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30195-2
