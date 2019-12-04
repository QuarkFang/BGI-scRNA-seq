import os
import re
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from multiprocessing import Process, Manager


def read_10x_mtx(filepath, multiprocess=True):
    dirs = os.listdir(filepath)
    list_dir = []
    for dir in dirs:
        if re.match(r'data(.*)', dir) is not None:
            list_dir.append(filepath + dir)
    genes, barcodes = read_tsv(list_dir)
    if multiprocess is True:
        matrix = read_mtx_multiprocess(list_dir)
        matrix = list(matrix)
    else:
        matrix = read_mtx_no_multiprocess(list_dir)
    return matrix, genes, barcodes


def read_tsv(list):
    genes = []
    barcodes = []
    for dir in list:
        genes.append(pd.read_csv(dir+'/genes.tsv', header=None, sep='\t'))
        barcodes.append(pd.read_csv(dir + '/barcodes.tsv', header=None, sep='\t'))
    return genes, barcodes


def read_mtx(list, matrix_list):
    for dir in list:
        filename = dir + '/matrix.mtx'
        mtx = mmread(filename).astype('float32')
        matrix_list.append(csr_matrix(mtx))


def read_mtx_multiprocess(list, process=2):
    manager = Manager()
    matrix_list = manager.list()
    jobs = []
    list_multiprocess = [[] for i in range(process)]
    for i, dir in enumerate(list):
        list_multiprocess[i%process].append(dir)

    for j in range(process):
        p = Process(target=read_mtx, args=(list_multiprocess[j], matrix_list))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()
    return matrix_list


def read_mtx_no_multiprocess(list):
    matrix_list = []
    for dir in list:
        filename = dir + '/matrix.mtx'
        mtx = mmread(filename).astype('float32')
        matrix_list.append(csr_matrix(mtx))
    return matrix_list


if __name__ == '__main__':
    # time_start = time()
    # result_multi = read_mtx_multiprocess(None)
    # time_end = time()
    # print(time_end-time_start)
    # time_start = time()
    # result_single = read_mtx_no_multiprocess(None)
    # time_end = time()
    # print(time_end - time_start)
    read_10x_mtx('./data/', multiprocess=True)

