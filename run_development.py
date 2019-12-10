import numpy as np
import pandas as pd
import scanpy as sc


if __name__ == '__main__':
    sc.settings.verbosity = 3
    sc.logging.print_versions()

    sc.settings.set_figure_params(dpi=80)

    adata = sc.read_h5ad('./cache/data.h5ad')

    adata.var_names_make_unique()

    # sc.pl.highest_expr_genes(adata, n_top=20)

    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=3)

    mito_genes = adata.var_names.str.startswith('MT-')

    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

    adata.obs['n_counts'] = adata.X.sum(axis=1).A1

    adata = adata[adata.obs['n_genes'] < 15000, :]
    adata = adata[adata.obs['percent_mito'] < 0.05, :]

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

    sc.pp.log1p(adata)

    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=40)

    adata = adata[:, adata.var['highly_variable']]

    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'], n_jobs=None)

    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver='arpack')

    sc.pl.pca_variance_ratio(adata, log=True)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=8)

    a = adata.uns['neighbors']['connectivities'].toarray()

    sc.tl.umap(adata)

    sc.settings.set_figure_params(dpi=800)

    sc.tl.louvain(adata)
    sc.pl.umap(adata, color=['louvain'])
