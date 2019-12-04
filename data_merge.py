from read import read_10x_mtx
from save import save_10x_h5
from scipy.sparse import hstack, csr_matrix


def merge_10x_mtx(matrix, genes, barcodes, axis='rows'):
    if axis not in ['rows', 'cols']:
        raise ValueError("the 'axis' should be 'rows' or 'cols'")
    if axis is 'rows':
        matrix_stack = matrix[0]
        barcodes_stack = barcodes[0]
        for i in range(1, len(matrix)):
            matrix_stack = hstack((matrix_stack, matrix[i]))
            barcodes_stack = barcodes_stack.append(barcodes[i], ignore_index=True)
        genes_stack = genes[0]

        return csr_matrix(matrix_stack), genes_stack, barcodes_stack

    if axis is 'cols':
        raise ValueError("It's not supported now")


if __name__ == '__main__':
    matrix, genes, barcodes = read_10x_mtx('./data/', multiprocess=True)
    matrix_stack, genes_stack, barcodes_stack = merge_10x_mtx(matrix, genes, barcodes)
    save_10x_h5(matrix_stack, genes_stack, barcodes_stack)
