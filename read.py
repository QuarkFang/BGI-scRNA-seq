from scipy.io import mmread
from scipy.sparse import csr_matrix
from multiprocessing import Process, Manager
from time import time


def read():
    pass


def read_mtx(filepath, matrix_list):
    filename = filepath + '/matrix.mtx'
    mtx = mmread(filename).astype('float32')
    matrix_list.append(csr_matrix(mtx))


def read_mtx_multithread(list, thread=8):
    manager = Manager()
    matrix_list = manager.list()
    jobs = []
    for i in range(thread):
        filepath = './data/data%d' % i
        p = Process(target=read_mtx, args=(filepath, matrix_list))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()
    return matrix_list


def read_mtx_no_multithread(list, thread=8):
    matrix_list = []
    for i in range(thread):
        filepath = './data/data%d' % i
        read_mtx(filepath, matrix_list)


if __name__ == '__main__':
    time_start = time()
    read_mtx_multithread(None)
    time_end = time()
    print(time_end-time_start)
    time_start = time()
    read_mtx_no_multithread(None)
    time_end = time()
    print(time_end - time_start)

