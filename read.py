import os
import re
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from multiprocessing import Process, Manager
from time import time


def read_10x_mtx(filepath, multiprocess=2):
    dirs = os.listdir(filepath)
    list_dir = []
    for dir in dirs:
        if re.match(r'data(.*)', dir) is not None:
            list_dir.append(filepath + dir)
    genes, barcodes = read_tsv(list_dir)
    if multiprocess > 1:
        matrix = read_mtx_multiprocess(list_dir, multiprocess)
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


def read_mtx(list, matrix_dict):
    for dir in list:
        index = dir[-1]
        filename = dir + '/matrix.mtx'
        mtx = mmread(filename).astype('float32')
        matrix_dict[index] = csr_matrix(mtx)


def read_mtx_multiprocess(dir_list, process=2):
    manager = Manager()
    matrix_dict = manager.dict()
    jobs = []
    list_multiprocess = [[] for i in range(process)]
    for i, dir in enumerate(dir_list):
        list_multiprocess[i%process].append(dir)

    for j in range(process):
        p = Process(target=read_mtx, args=(list_multiprocess[j], matrix_dict))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()

    keys = matrix_dict.keys()
    keys.sort()
    matrix_list = [matrix_dict[key] for key in keys]

    return matrix_list


def read_mtx_no_multiprocess(list):
    matrix_list = []
    for dir in list:
        filename = dir + '/matrix.mtx'
        mtx = mmread(filename).astype('float32')
        matrix_list.append(csr_matrix(mtx))
    return matrix_list


if __name__ == '__main__':
    time_start = time()
    read_10x_mtx('./data/', multiprocess=10)
    time_end = time()
    print(time_end-time_start)
    time_start = time()
    read_10x_mtx('./data/', multiprocess=1)
    time_end = time()
    print(time_end - time_start)
    # read_10x_mtx('./data/', multiprocess=2)

