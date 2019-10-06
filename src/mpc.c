#include "mpc.h"
#include "block.h"
#include "lin_alg.h"

// TODO: Imlement const correctness!

/* A convenience function, the implementation is
 * useful only in this context, as it does not check
 * for zero elements.
 * Takes ownership of x.
 */
csc* dns_to_csc(c_int m, c_int n, c_float *x)
{
  csc *M = csc_spalloc(m, n, m*n, 0, 0);

  M->p[0] = 0;
  int i = 0, j = 0;
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      M->i[j*m+i] = i;
      M->p[j+1] = M->p[j]+m; 
    }
  }
  M->x = x;
}

csc* mpc_to_osqp_P(const csc* Q,
                   const csc* QN,
                   const csc* R, c_int N)
{
  csc* P = OSQP_NULL;

  csc* I = csc_eye(N, N, 0);
  csc* P_Q = csc_kron(I, Q);
  csc* P_R = csc_kron(I, R);
  P = csc_block_diag(3, P_Q, QN, P_R);

  csc_spfree(I);
  csc_spfree(P_Q);
  csc_spfree(P_R);
  
  return P;
}

c_float* mpc_to_osqp_q_const_xr(const csc* Q,
                                const csc* QN,
                                const c_float* xr,
                                c_int nu, c_int N)
{
  c_float *q = OSQP_NULL;
  c_float *y1 = c_malloc(Q->m * sizeof(c_float));
  c_float *y2 = c_malloc(QN->m * sizeof(c_float));
  
  mat_vec(Q, xr, y1, 0);
  csc *minus_Q_dot_xr = dns_to_csc(Q->m, 1, y1);
  mat_mult_scalar(minus_Q_dot_xr, -1);
  csc *O = csc_ones(N, 1);
  csc *qQ = csc_kron(O, minus_Q_dot_xr);
  csc_spfree(O);
  csc_spfree(minus_Q_dot_xr);
  
  mat_vec(QN, xr, y2, 0);
  csc *minus_QN_dot_xr = dns_to_csc(QN->m, 1, y2);
  mat_mult_scalar(minus_QN_dot_xr, -1);
  csc *Z = csc_zeros(N*nu, 1);

  csc *q_csc = csc_vstack(3, qQ, minus_QN_dot_xr, Z);
  q = csc_to_dns(q_csc);

  csc_spfree(q_csc);
  csc_spfree(Z);
  csc_spfree(minus_QN_dot_xr);
  csc_spfree(qQ);
  
  return q;
}

c_float* mpc_to_osqp_q(csc* Q,
                       csc* QN,
                       const c_float* xr,
                       c_int nu, c_int N)
{
  c_int nx = Q->m;
  c_float *q = c_calloc((N+1) * nx + N * nu, sizeof(c_float));

  csc *I = csc_eye(N, N, 0);
  mat_mult_scalar(Q, -1.0);
  csc* repQ = csc_kron(I, Q);

  // (-1)*(-Q) = Q
  // We do this immediately, because the user
  // may pass Q as QN
  mat_mult_scalar(Q, -1.0);
  csc_spfree(I);
  
  mat_vec(repQ, xr, q, 0);

  csc_spfree(repQ);

  mat_mult_scalar(QN, -1.0);
  mat_vec(QN, &xr[N*nx], &q[N*nx], 0);

  mat_mult_scalar(QN, -1.0);
  
  return q;
}

csc* mpc_to_osqp_A(const csc* Ad,
                   const csc* Bd,
                   c_int N)
{
  csc* A = OSQP_NULL;
  c_int nx = Ad->n;
  c_int nu = Bd->n;

  csc* IN = csc_eye(N+1, N+1, 0);
  csc* minus_Inx = csc_eye(nx, nx, 0);
  mat_mult_scalar(minus_Inx, -1);
  csc* Ax_I = csc_kron(IN, minus_Inx);
  csc_spfree(IN);
  csc_spfree(minus_Inx);

  csc* IN_k1 = csc_eye(N+1, N+1, -1);
  csc* Ax_Ad = csc_kron(IN_k1, Ad);
  csc_spfree(IN_k1);

  csc* Ax = csc_add(Ax_I, Ax_Ad, 1, 1);

  csc_spfree(Ax_I);
  csc_spfree(Ax_Ad);
  
  csc* Z = csc_zeros(1,N);
  IN = csc_eye(N, N, 0);
  csc* Z_IN = csc_vstack(2, Z, IN);
  csc* Bu = csc_kron(Z_IN, Bd);
 
  csc_spfree(Z);
  csc_spfree(IN);
  csc_spfree(Z_IN);

  csc *Aeq = csc_hstack(2, Ax, Bu);

  csc_spfree(Ax);
  csc_spfree(Bu);

  csc *Aineq = csc_eye((N+1)*nx + N*nu, (N+1)*nx + N*nu, 0);
  A = csc_vstack(2, Aeq, Aineq);

  csc_spfree(Aeq);
  csc_spfree(Aineq);
  
  return A;
}

csc* lpv_a_mpc_to_osqp_A(csc *Ad[], const csc* Bd, const csc *Ineq[], c_int N)
{
  csc* A = OSQP_NULL;
  c_int nx = Ad[0]->n;
  c_int nu = Bd->n;
  c_int k = 0, j = 0;

  csc* IN = csc_eye(N+1, N+1, 0);
  csc* minus_Inx = csc_eye(nx, nx, 0);
  mat_mult_scalar(minus_Inx, -1);
  csc* Ax_I = csc_kron(IN, minus_Inx);
  csc_spfree(IN);
  csc_spfree(minus_Inx);

  c_int nzmax = Ad[0]->nzmax;
  csc* Ax_Ad = csc_spalloc((N+1)*nx, (N+1)*nx, N*nzmax, 1, 0);
  for (k = 0; k < N; k++)
  {
    for (j = 0; j < nzmax; j++)
    {
      Ax_Ad->x[k*nzmax+j] = Ad[k]->x[j];
      Ax_Ad->i[k*nzmax+j] = Ad[k]->i[j] + (k+1)*nx;
    }
    for (j = 0; j <= nx; j++)
    {
      Ax_Ad->p[k*nx+j] = Ad[k]->p[j] + k*nzmax;
    }
  }
  // Last nx columns of Ax_Ad are zeros
  for (j = N*nx; j <= (N+1)*nx; j++)
  {
    Ax_Ad->p[j] = N*nzmax;
  }

  csc* Ax = csc_add(Ax_I, Ax_Ad, 1, 1);

  csc_spfree(Ax_I);
  csc_spfree(Ax_Ad);

  csc* Z = csc_zeros(1,N);
  IN = csc_eye(N, N, 0);
  csc* Z_IN = csc_vstack(2, Z, IN);
  csc* Bu = csc_kron(Z_IN, Bd);

  csc_spfree(Z);
  csc_spfree(IN);
  csc_spfree(Z_IN);

  csc *Aeq = csc_hstack(2, Ax, Bu);

  csc_spfree(Ax);
  csc_spfree(Bu);

  csc *Aineq = OSQP_NULL;
  if (Ineq == OSQP_NULL)
  {
    // No additional input constraints.
    Aineq = csc_eye((N+1)*nx + N*nu, (N+1)*nx + N*nu, 0);
  }
  else
  {
    // TODO: Include Ineq into Aineq!!!
    Aineq = csc_eye((N+1)*nx + N*nu, (N+1)*nx + N*nu, 0);
  }
  A = csc_vstack(2, Aeq, Aineq);

  csc_spfree(Aeq);
  csc_spfree(Aineq);

  return A;
}

c_float* mpc_to_osqp_bound(const c_float* minus_x0,
                           const c_float* x_bound,
                           const c_float* xN_bound,
                           const c_float* u_bound,
                           const c_float* other_bound,
                           c_int nx,
                           c_int nu,
                           c_int N,
                           c_int n_other)
{
  c_float* vec_z = vec_zeros(N*nx);
  c_float* beq = vec_cat(nx, minus_x0, N*nx, vec_z);
  c_free(vec_z);

  c_float* bineq_x = vec_rep(nx, x_bound, N);
  c_float* bineq_x_xN = vec_cat(N*nx, bineq_x, nx, xN_bound);
  c_float* bineq_u = vec_rep(nu, u_bound, N);
  c_float* bineq = vec_cat((N+1)*nx, bineq_x_xN, N*nu, bineq_u);

  c_free(bineq_u);
  c_free(bineq_x_xN);
  c_free(bineq_x);

  c_float* bound = vec_cat((N+1)*nx, beq, (N+1)*nx+N*nu, bineq);

  c_free(beq);
  c_free(bineq);

  c_float* bound_full = bound;
  if (other_bound != OSQP_NULL)
  {
    bound_full = vec_cat(2*(N+1)*nx+N*nu, bound, n_other, other_bound);
    c_free(bound);
  }

  return bound_full;
}

