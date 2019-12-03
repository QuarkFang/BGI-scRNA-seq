import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from scipy.io import mmwrite
from scipy.sparse import csr_matrix


def read_tsv(filepath):
    return pd.read_csv(filepath, header=None, sep='\t')


def generation(filepath, filenum, rows=None, cols=None):
    if rows is None and cols is None:
        raise ValueError("at least one of 'rows' and 'cols' needs value")
    elif rows is None and cols is not None:
        rows = len(read_tsv(filepath + 'genes.tsv'))
        if not isinstance(cols, int):
            raise ValueError("'cols' should be int")
        data = generation_data_rows(filepath, filenum, rows, cols)
        write_data(data)
    elif rows is not None and cols is None:
        raise ValueError("It's not supported now")
    else:
        raise ValueError("It's not supported now")


def generation_data_rows(filepath, filenum, rows, cols):
    data = []
    genes = read_tsv(filepath + '/genes.tsv')
    for i in range(filenum):
        matrix = generation_mtx(rows, cols)
        barcodes = read_tsv(filepath + '/barcodes.tsv')[i*cols:(i+1)*cols]
        data.append((matrix, genes, barcodes))
    return data


def generation_mtx(rows, cols):
    return csr_matrix(np.random.randint(0, 2, size=[rows, cols]))


def write_data(data):
    for single in tqdm(enumerate(data)):
        (i, (matrix, genes, barcodes)) = single
        dir = './data/data'+str(i)
        os.mkdir(dir)
        mmwrite(dir + '/matrix.mtx', matrix)
        genes.to_csv(dir + '/genes.tsv', sep='\t', index=False, header=False)
        barcodes.to_csv(dir + '/barcodes.tsv', sep='\t', index=False, header=False)


if __name__ == '__main__':
    generation('./data/tsv/', 10, cols=2000)
