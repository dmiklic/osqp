/* Operations on block matrices */

#ifndef BLOCK_H
#define BLOCK_H

#include "types.h"
#include "cs.h"

#include <stdarg.h>

/**
 * Create Compressed-Column-Sparse zero matrix
 * Calls csc_spalloc to allocate memory, the caller is 
 * responsible for freeing memory with csc_spfree.
 * @param  m     Number of rows
 * @param  n     Number of columns
 * @return       New matrix pointer
 */
csc* csc_zeros(c_int m, c_int n);

/**
 * Create Compressed-Column-Sparse matrix of ones
 * Calls csc_spalloc to allocate memory, the caller is 
 * responsible for freeing memory with csc_spfree.
 * @param  m     Number of rows
 * @param  n     Number of columns
 * @return       New matrix pointer
 */
csc* csc_ones(c_int m, c_int n);

/**
 * Create Compressed-Column-Sparse square identity matrix
 * Calls csc_spalloc to allocate memory, the caller is 
 * responsible for freeing memory with csc_spfree.
 * @param  n     Matrix size
 * @return       New matrix pointer
 */
csc* csc_eye(c_int m, c_int n, c_int k);

/**
 * Create Compressed-Column-Sparse square diagonal matrix
 * Calls csc_spalloc to allocate memory, the caller is 
 * responsible for freeing memory with csc_spfree.
 * Makes a copy of elems.
 * @param  n     Matrix size
 * @param  elems Elements on the diagonal
 * @return       New matrix pointer
 */
csc* csc_diag(c_int m, c_int n, const c_float *elems, c_int k);

/**
 * Transpose a Compressed-Column-Sparse matrix.
 * @param  M Matrix to transpose.
 * @return Returns a new matrix T = M'
 *         T is allocated using csc_spalloc and must be
 *         freed using csc_spfree
 */
csc* csc_tpose(const csc* M);

/**
 * Stack Compressed-Column-Sparse matrices horizontally.
 * Calls csc_spalloc to allocate memory, the caller is 
 * responsible for freeing memory with csc_spfree.
 * @param  count Number of matrices to stack.
 * @param  List of matrices to stack.
 * @return Pointer to the newly crated matrix or OSQP_NULL if 
 *         matrix dimensions don't match or there is insufficient
 *         memory.
 */
csc* csc_hstack_list(c_int count, csc *Mlist[]);

/**
 * Stack Compressed-Column-Sparse matrices horizontally.
 * Calls csc_spalloc to allocate memory, the caller is 
 * responsible for freeing memory with csc_spfree.
 * @param  count Number of matrices to stack.
 * @param  Matrices to stack (all optional arguments must be csc*)
 * @return Pointer to the newly crated matrix or OSQP_NULL if 
 *         matrix dimensions don't match or there is insufficient
 *         memory.
 */
csc* csc_hstack(c_int count, ...);

/**
 * Stack Compressed-Column-Sparse matrices vertically.
 * Calls csc_spalloc to allocate memory, the caller is 
 * responsible for freeing memory with csc_spfree.
 * @param  count Number of matrices to stack.
 * @param  Matrices to stack (all optional arguments must be csc*)
 * @return Pointer to the newly crated matrix.
 */
csc* csc_vstack(c_int count, ...);

/**
 * Stack Compressed-Column-Sparse matrices vertically.
 * Calls csc_spalloc to allocate memory, the caller is
 * responsible for freeing memory with csc_spfree.
 * @param  count Number of matrices to stack.
 * @param  List of matrices to stack.
 * @return Pointer to the newly crated matrix.
 */
csc* csc_block_diag_list(c_int count, csc *Mlist[]);

/**
 * Stack Compressed-Column-Sparse matrices vertically.
 * Calls csc_spalloc to allocate memory, the caller is
 * responsible for freeing memory with csc_spfree.
 * @param  count Number of matrices to stack.
 * @param  Matrices to stack (all optional arguments must be csc*)
 * @return Pointer to the newly crated matrix.
 */
csc* csc_block_diag(c_int count, ...);

/**
 * Kronecker product of Compressed-Column-Sparse matrices.
 * Calls csc_spalloc to allocate memory, the caller is
 * responsible for freeing memory with csc_spfree.
 * @param  A Pointer to the first operand
 * @param  B Pointer to the second operand
 * @return Pointer to the newly crated matrix.
 */
csc* csc_kron(const csc *A, const csc *B);

//void csc_set(c_int i, c_int j, c_float val);
//c_float csc_get(c_int i, c_int j);
void print_csc_matrix_as_dns(csc* M, const char *name);

/*
 * Add Compressed-Column-Sparse matrices.
 * Computes C = alpha * A + beta * B
 * Allocates memory for the result. The caller is
 * responsible for freeing the memory by calling
 * csc_spfree.
 * Does not perform "compression", i.e. if the addition
 * produces a zero element, it is not removed from the 
 * sparse representation(?).
 * @param  A Pointer to the first operand
 * @param  B Pointer to the second operand
 * @param  alpha
 * @param  beta
 * @return The result C = alpha * A + beta * B
 */
csc* csc_add(const csc *A, const csc *B, double alpha, double beta);
#endif
