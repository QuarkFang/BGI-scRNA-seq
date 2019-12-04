from anndata import AnnData
import os


def save_10x_h5(matrix_stack, genes_stack, barcodes_stack, dtype: str='float32'):
    adata = AnnData(matrix_stack, dtype=dtype).T
    var_names = genes_stack[1]
    adata.var_names = var_names
    adata.var['gene_ids'] = genes_stack[0].values
    adata.obs_names = barcodes_stack[0]

    path = os.fspath('./cache/data.h5ad')
    adata.write(path)
