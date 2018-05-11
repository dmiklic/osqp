/* Implementation of block operations on sparse matrices
 *
 * NOTE: The implementations are not optimized in terms of performance 
 * or memory consumption. Their purpose is to simplify QP formulation
 * of MPC problems.
 *
 */

#include "block.h"
#include "cs.h"

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
  /*
  c_int i = 0;
  csc *I = csc_spalloc(n, n, n, 1, 0);
  for (i = 0; i < n; i++)
  {
    I->x[i] = 1.00000000000000000000;
    I->i[i] = i;
    I->p[i] = i;
  }
  I->p[n] = n;
  */
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

//void csc_set(c_int i, c_int j, c_float val);
//void csc_get(c_int i, c_int j);

void print_csc_matrix_as_dns(csc* M, const char* name)
{
  c_float* Mdns = csc_to_dns(M);
  print_dns_matrix(Mdns, M->m, M->n, name);
  c_free(Mdns);
}
