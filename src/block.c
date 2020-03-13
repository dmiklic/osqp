/* Implementation of block operations on sparse matrices
 *
 * NOTE: The implementations are not optimized in terms of performance 
 * or memory consumption. Their purpose is to simplify QP formulation
 * of MPC problems.
 *
 */

#include "block.h"
#include "cs.h"
#include "util.h"

#include <limits.h>

csc* csc_zeros(c_int m, c_int n)
{
  csc *Z = c_calloc(1, sizeof(csc));
  Z->m = m;
  Z->n = n;
  Z->nzmax = 0;
  Z->nz = -1;
  Z->p = c_calloc(n+1, sizeof(c_int));
  Z->i = OSQP_NULL;
  Z->x = OSQP_NULL;
  return Z;
}

csc* csc_ones(c_int m, c_int n)
{
  csc *O = csc_spalloc(m, n, m*n, 1, 0);
  O->p[0] = 0;
  int i = 0, j = 0;
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      O->x[j*m+i] = 1.00000000000000000000;
      O->i[j*m+i] = i;
      O->p[j+1] = O->p[j]+m; 
    }
  }
  return O;
}

csc* csc_eye(c_int m, c_int n, c_int k)
{
  c_int nnz = 0, i = 0, j = 0;
  csc *I = OSQP_NULL;
  if (k >= 0)
  {
    nnz = c_min(m, n-k);
    I = csc_spalloc(m, n, nnz, 1, 0);
    for (i = 0; i < nnz; i++)
    {
      I->i[i] = i;
      I->x[i] = 1.0;
    }
    for (j = 0; j <= k; j++)
    {
      I->p[j] = 0;
    }
    for (j = k+1; j < n+1; j++)
    {
      I->p[j] = c_min(j - k, nnz);
    }
  }
  else if (k < 0)
  {
    nnz = c_min(m+k, n);
    I = csc_spalloc(m, n, nnz, 1, 0);
    for (i = 0; i < nnz; i++)
    {
      I->i[i] = i - k;
      I->x[i] = 1.0;
    }
    for (j = 0; j < n+1; j++)
    {
      I->p[j] = c_min(j, nnz);
    }
  }
  return I;
}

csc* csc_diag(c_int m, c_int n, const c_float *elems, c_int k)
{
  c_int i = 0;
  csc *D = csc_eye(m, n, k);
  for (i = 0; i < D->nzmax; i++)
  {
    D->x[i] = elems[i];
  }
  return D;
}

// Alorithm described here:
// https://stackoverflow.com/questions/16721863/how-does-matlab-transpose-a-sparse-matrix
// TODO: Remove temp
csc* csc_tpose(const csc* M)
{
  // Special handling of zero matrices
  if (M->nzmax == 0)
  {
    return csc_zeros(M->n, M->m);
  }
  csc* T = csc_spalloc(M->n, M->m, M->nzmax, 1, 0);
  c_int* temp = c_calloc(M->m+1, sizeof(c_int));
  int i, j;
  // Count number of elements in each row of A
  // (which are going to become columns of T)
  for (i = 0; i < M->n; ++i)
  {
    for (j = M->p[i]; j < M->p[i+1]; ++j)
    {
      ++temp[M->i[j]+1];
    }
  }
  // T->p is the cumulative sum of the number
  // of elements in each column of T
  T->p[0] = 0;
  for (i = 1; i <= M->m; ++i)
  {
    T->p[i] = T->p[i-1] + temp[i];
    temp[i] = T->p[i];
  }
  for (i = 0; i < M->n; ++i)
  {
    for (j = M->p[i]; j < M->p[i+1]; ++j)
    {
      T->i[temp[M->i[j]]] = i;
      T->x[temp[M->i[j]]++] = M->x[j];
    }
  }
  c_free(temp);
  return T;
}

csc* csc_hstack_list(c_int count, csc *Mlist[])
{
  csc *H = OSQP_NULL;
  c_int i = 0;
  c_int numrows_ok = 1;
  
  int m = Mlist[0]->m;
  int n = Mlist[0]->n;
  int nzmax = Mlist[0]->nzmax;
  // Check if all row numbers match
  for (i = 1; i < count; i++)
  {
    if (Mlist[i]->m == m)
    {
      n += Mlist[i]->n;
      nzmax += Mlist[i]->nzmax;
    }
    else
    {
      numrows_ok = 0;
      break;
    }
  }

  if (numrows_ok)
  {
    H = csc_spalloc(m, n, nzmax, 1, 0); 
    // Reset variables to keep track of the
    // number of elements copied into H
    
    nzmax = 0;
    n = 0;
    for (i = 0; i < count; i++)
    {
      int k;
      for (k = 0; k < Mlist[i]->nzmax; k++)
      {
        H->x[nzmax+k] = Mlist[i]->x[k];
        H->i[nzmax+k] = Mlist[i]->i[k];
      }
      // TODO: Replace with memcpy(?)
      for (k = 0; k < Mlist[i]->n; k++)
      {
        H->p[n+k] = nzmax + Mlist[i]->p[k];
      }
      nzmax += Mlist[i]->nzmax;
      n += Mlist[i]->n;
    }
    H->p[n] = nzmax;
  }
  return H;
}

csc* csc_hstack(c_int count, ...)
{
  csc **Mlist = c_calloc(count, sizeof(csc*));
  va_list args;
  va_start(args, count);
  c_int i = 0;
  for (i = 0; i < count; i++)
  {
    Mlist[i] = va_arg(args, csc*);
  }
  va_end(args);
  csc *H = csc_hstack_list(count, Mlist);
  c_free(Mlist);
  return H;
}

csc* csc_vstack_list(c_int count, csc *Mlist[])
{
  // We transpose everything, stack horizontally
  // and then transpose back.
  
  csc *V = OSQP_NULL;
  
  // Storage for transposed matrices
  csc **Tlist = c_calloc(count, sizeof(csc*));

  int i;
  for (i = 0; i < count; i++)
  {
    Tlist[i] = csc_tpose(Mlist[i]);
  }

  csc *H = csc_hstack_list(count, Tlist);
  V = csc_tpose(H);

  csc_spfree(H);
  for (i = 0; i < count; i++)
  {
    csc_spfree(Tlist[i]);
  }
  c_free(Tlist);
  
  return V;
}

csc* csc_vstack(c_int count, ...)
{   
  // Storage for transposed matrices
  csc **Mlist = c_calloc(count, sizeof(csc*));

  va_list args;
  va_start(args, count);
  int i;
  for (i = 0; i < count; i++)
  {
    Mlist[i] = va_arg(args, csc*);
  }
  va_end(args);

  csc *V = csc_vstack_list(count, Mlist);
  c_free(Mlist);
  
  return V;
}

csc* csc_block_diag_list(c_int count, csc *Mlist[])
{
  csc *B = OSQP_NULL;
  c_int *n_cum = c_calloc(count+2, sizeof(c_int));
  if (count == 1)
  {
    B = copy_csc_mat(Mlist[0]);
  }
  else if (count > 1)
  {
    int i = 0;
    for (i = 0; i < count; i++)
    {
      n_cum[i+1] = n_cum[i] + Mlist[i]->n;
    }
    n_cum[count+1] = n_cum[count];
    csc **Brows = c_calloc(count, sizeof(csc*));
    for (i = 0; i < count; i++)
    {
      csc* Z1 = csc_zeros(Mlist[i]->m, n_cum[i]);
      csc* Z2 = csc_zeros(Mlist[i]->m, n_cum[count+1]-n_cum[i+1]);
      Brows[i] = csc_hstack(3, Z1, Mlist[i], Z2);
      csc_spfree(Z1);
      csc_spfree(Z2);
    }
    B =  csc_vstack_list(count, Brows);

    for (i = 0; i < count; i++)
    {
      csc_spfree(Brows[i]);
    }
    c_free(Brows);
    c_free(n_cum);
  }
  return B;
}

csc* csc_block_diag(c_int count, ...)
{
  csc **Mlist = c_calloc(count, sizeof(csc*));
  va_list args;
  va_start(args, count);
  c_int i = 0;
  for (i = 0; i < count; i++)
  {
    Mlist[i] = va_arg(args, csc*);
  }
  va_end(args);
  csc *B = csc_block_diag_list(count, Mlist);
  c_free(Mlist);
  return B;
}

csc* csc_kron(const csc *A, const csc *B)
{
  csc* K = csc_spalloc(A->m * B->m,A->n * B->n, A->nzmax * B->nzmax, 1, 0);

  K->p[0] = 0;
  c_int idxK = 0;
  c_int jA = 0, jB = 0;
  for (jA = 0; jA < A->n; jA++)
  {
    for (jB = 0; jB < B->n; jB++)
    {
      K->p[jA*B->n+jB+1] = K->p[jA*B->n+jB]
        + (A->p[jA+1]-A->p[jA]) * (B->p[jB+1]-B->p[jB]);
      c_int idxA = 0, idxB = 0;
      for (idxA = A->p[jA]; idxA < A->p[jA+1]; idxA++)
      {
        for (idxB = B->p[jB]; idxB < B->p[jB+1]; idxB++)
        {
          K->x[idxK] = A->x[idxA] * B->x[idxB];
          K->i[idxK++] = A->i[idxA] * B->m + B->i[idxB];
        }
      }
    }
  }
  return K;
}

void print_csc_matrix_as_dns(csc* M, const char* name)
{
  c_float* Mdns = csc_to_dns(M);
  print_dns_matrix(Mdns, M->m, M->n, name);
  c_free(Mdns);
}

c_float* vec_zeros(c_int n)
{
  return c_calloc(n, sizeof(c_float));
}

c_float* vec_ones(c_int n)
{
  c_float* vec = c_malloc(n * sizeof(c_float));
  c_int i = 0;
  for (i = 0; i < n; i++)
  {
    vec[i] = 1.0;
  }
  return vec;
}

c_float* vec_cat(c_int n1, const c_float* v1, c_int n2, const c_float* v2)
{
  c_float* vec = c_malloc((n1+n2) * sizeof(c_float));
  memcpy(vec, v1, n1*sizeof(c_float));
  memcpy(vec+n1, v2, n2*sizeof(c_float));
  return vec;
}

c_float* vec_rep(c_int n, const c_float* v, c_int count)
{
  c_float* vec = c_malloc(n*count*sizeof(c_float));
  int i = 0;
  for (i = 0; i < count; i++)
  {
    memcpy(vec+i*n, v, n*sizeof(c_float));
  }
  return vec;
}

//void csc_set(c_int i, c_int j, c_float val);
//void csc_get(c_int i, c_int j);

// All code below taken directly from:
// http://people.sc.fsu.edu/~jburkardt/c_src/csparse/csparse.html

#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))

/* free workspace and return a sparse matrix result */
csc *cs_done (csc *C, void *w, void *x, c_int ok)
{
    c_free (w) ;			/* free workspace */
    c_free (x) ;
    return (ok ? C : csc_spfree (C)) ;	/* return result if OK, else free it */
}

/* wrapper for realloc */
void *cs_realloc (void *p, c_int n, size_t size, c_int *ok)
{
    void *p2 ;
    *ok = !CS_OVERFLOW (n,size) ;	    /* guard against int overflow */
    if (!(*ok)) return (p) ;		    /* p unchanged if n too large */
    p2 = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
    *ok = (p2 != NULL) ;
    return ((*ok) ? p2 : p) ;		    /* return original p if failure */
}

/* change the max # of entries sparse matrix */
c_int cs_sprealloc (csc *A, c_int nzmax)
{
    c_int ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    nzmax = (nzmax <= 0) ? (A->p [A->n]) : nzmax ;
    A->i = cs_realloc (A->i, nzmax, sizeof (c_int), &oki) ;
    if (A->nz >= 0) A->p = cs_realloc (A->p, nzmax, sizeof (c_int), &okj) ;
    if (A->x) A->x = cs_realloc (A->x, nzmax, sizeof (c_float), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
c_int cs_scatter (const csc *A, c_int j, c_float beta, c_int *w, c_float *x, c_int mark,
                  csc *C, c_int nz)
{
    c_int i, p, *Ap, *Ai, *Ci ;
    c_float *Ax ;
    if (!A || !w || !C) return (-1) ;		/* ensure inputs are valid */
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
	i = Ai [p] ;				/* A(i,j) is nonzero */
	if (w [i] < mark)
	{
	    w [i] = mark ;			/* i is new entry in column j */
	    Ci [nz++] = i ;			/* add i to pattern of C(:,j) */
	    if (x) x [i] = beta * Ax [p] ;	/* x(i) = beta*A(i,j) */
	}
	else if (x) x [i] += beta * Ax [p] ;	/* i exists in C(:,j) already */
    }
    return (nz) ;
}

csc *csc_add ( const csc *A, const csc *B, double alpha, double beta )
/*
  Purpose:

    CS_ADD computes C = alpha*A + beta*B for sparse A and B.

  Reference:

    Timothy Davis,
    Direct Methods for Sparse Linear Systems,
    SIAM, Philadelphia, 2006.
*/
{
    c_int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
    c_float *x, *Bx, *Cx ;
    csc *C ;
    if (!A || !B) return (NULL) ;	/* check inputs */
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
    w = c_calloc (m, sizeof (c_int)) ;
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? c_malloc (m * sizeof (c_float)) : NULL ;
    C = csc_spalloc (m, n, anz + bnz, values, 0) ;
    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)
    {
	Cp [j] = nz ;			/* column j of C starts here */
	nz = cs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
	nz = cs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
	if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;			/* finalize the last column of C */
    cs_sprealloc (C, 0) ;		/* remove extra space from C */
    return (cs_done (C, w, x, 1)) ;	/* success; free workspace, return C */
}

