# Unsupervised Learning

The unsupervised part of the project focuses on understanding the intrinsic structure of the data without using target labels. For each dataset (`unsupervised_smartseq.ipynb`, `unsupervised_dropseq.ipynb`) the same pipeline runs independently on two cell lines (**MCF7** and **HCC1806**), and the recovered structure is compared back to the held-out Hypo/Norm condition to decide whether hypoxia is the dominant transcriptional axis.

## Dimensionality reduction

The following methods were used for low-dimensional embedding and visualization, within a `scanpy` pipeline:

- **PCA**
- **UMAP**

PCA (30 components) builds the latent space; a nearest-neighbour graph on the PCA space feeds **UMAP** for 2D visualization. Embeddings were explored to assess how clearly different cell populations separate in latent space — coloured first by the held-out condition, then by the recovered clustering.

A key evaluation metric used in this stage is the **silhouette score**, which measures how well-separated the clusters are in the embedding. We also locate where the signal lives by correlating each principal component with the binary condition (point-biserial *r*): on MCF7 hypoxia loads heavily on **PC1** (|r| up to 0.93 on SmartSeq), confirming it is the dominant axis, whereas on HCC1806 it leaks onto minor PCs.

## Clustering methods

The following clustering approaches were explored:

- **K-Means**
- **Leiden**

### K-Means
K-Means was applied on the PCA space, with the number of clusters `k` chosen by a silhouette sweep. Performance was analysed under different `k` values. A key caveat surfaced: the silhouette-optimal partition is not always biologically meaningful — on HCC1806 silhouette strongly prefers a `k=2` split that turns out to be unrelated to the condition (ARI ≈ 0), and on MCF7-DropSeq it prefers `k=3`, which over-splits the two real groups.

### Leiden
Leiden community detection was applied on the neighbour graph, with the **resolution** chosen by a silhouette sweep. Silhouette peaks at the lowest resolutions (2 clusters) for both lines; higher resolutions fragment the data and the score falls. Leiden proved the more reliable method here, recovering the condition well on MCF7 (ARI ≈ 0.98 SmartSeq, ≈ 0.88 DropSeq) where K-Means over-split.

Recovered clusterings were scored against the true labels with **Adjusted Rand Index (ARI)** and **Normalized Mutual Information (NMI)**, supported by contingency tables. Two further checks were run: a **hypoxia signature score** (`sc.tl.score_genes` on canonical hypoxia genes) as a biological sanity check, and a **`log1p` comparison** to test whether the missing log-transform limited recovery.

## Main takeaway from unsupervised analysis

The unsupervised experiments showed that:

- low-dimensional embeddings, especially **UMAP**, often provide visually meaningful separation between groups — clearly so for **MCF7**, where hypoxia is the dominant axis of variation and label-free clustering essentially reconstructs the experiment
- clustering performance depends heavily on preprocessing choices (normalisation, `log1p`) and data sparsity
- silhouette-optimal partitions are not always biologically meaningful, so ground-truth comparison (ARI/NMI) is essential
- the hypoxia signal is genuinely present even where geometric clustering fails (**HCC1806**, confirmed by the signature score), but it shares the variance structure with other factors
- **SmartSeq** outperforms **DropSeq** throughout: its deeper, full-length reads give a cleaner signal than DropSeq's shallow, dropout-heavy counts, where even 30 PCs capture only ~5 % of the variance
