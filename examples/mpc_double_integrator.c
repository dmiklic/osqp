#include "stdio.h"
#include "osqp.h"

#include "lin_alg.h"
#include "block.h"
#include "mpc.h"



int main(int argc, char **argv) {

  c_int nx = 2;
  c_int nu = 1;

  // System model

  c_int Ad_nnz = 3;
  c_float Ad_x[3] = { 1.0, 1.0, 1.0 };
  c_int Ad_i[3] = { 0, 0, 1 };
  c_int Ad_p[3] = { 0, 1, 3 };
  csc *Ad = csc_matrix(nx, nx, Ad_nnz, Ad_x, Ad_i, Ad_p);

  c_int Bd_nnz = 1;
  c_float Bd_x[1] = { 1.0 };
  c_int Bd_i[1] = { 1 };
  c_int Bd_p[2] = { 0, 1 };
  csc *Bd = csc_matrix(nx, nu, Bd_nnz, Bd_x, Bd_i, Bd_p);

  // Initial conditions and setpoints
  
  c_float x0[2] = { 3.0, -2.0 };
  c_float minus_x0[2] = { -3.0, 2.0 };
  c_float xr[2] = { 0.0, 0.0 };

  // Costs

  csc *Q = csc_eye(2,2,0);
  c_float QN_x[4] = { 2.668689083642158, 1.726606170754336,
                      1.726606170754336, 2.881168868886953 };
  c_int QN_nnz = 4;
  c_int QN_i[4] = { 0, 1, 0, 1 };
  c_int QN_p[3] = { 0, 2, 4 };
  csc *QN = csc_matrix(2, 2, QN_nnz, QN_x, QN_i, QN_p);
  csc *R = csc_eye(1, 1, 0);
  mat_mult_scalar(R, 0.1);

  // Constraints
  c_float x_min[2] = { -10.0, -10.0 };
  c_float xN_min[2] = { -HUGE_VAL, -HUGE_VAL };
  c_float u_min[1] = { -1.0 };
  c_float x_max[2] = { 10.0, 10.0 };
  c_float xN_max[2] = { HUGE_VAL, HUGE_VAL };
  c_float u_max[1] = { 1.0 };
  
  // Prediction horizon
  c_int N = 6;

  // Cast as QP problem
  csc *P = mpc_to_osqp_P(Q, QN, R, N);
  //print_csc_matrix_as_dns(P, "P");
  c_float *q = mpc_to_osqp_q(Q, QN, xr, nu, N);
  //print_dns_matrix(q, N*nx + nx + N*nu, 1, "q");
  csc *A = mpc_to_osqp_A(Ad, Bd, N);
  //print_csc_matrix_as_dns(A, "A");
  c_float *l = mpc_to_osqp_bound(minus_x0, x_min, xN_min, u_min, nx, nu, N);
  //print_dns_matrix(l, 1, (N+1)*nx + (N+1)*nx + N*nu, "l");
  c_float *u = mpc_to_osqp_bound(minus_x0, x_max, xN_max, u_max, nx, nu, N);
  //print_dns_matrix(u, 1, (N+1)*nx + (N+1)*nx + N*nu, "u");

  // Invoke OSQP
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  OSQPWorkspace *work;
  OSQPData *data;

  data = (OSQPData *)c_malloc(sizeof(OSQPData));
  data->n = (N+1)*nx + N*nu; // Number of variables
  data->m = (N+1)*nx + (N+1)*nx + N*nu; // Number of constraints
  data->P = P;
  data->q = q;
  data->A = A;
  data->l = l;
  data->u = u;

  osqp_set_default_settings(settings);
  work = osqp_setup(data, settings);

  osqp_solve(work);
  print_dns_matrix(work->solution->x + (N+1)*nx, 1, N*nu, "u*");
  // Clean up
  osqp_cleanup(work);
  
  // OSQP variables
  c_free(u);
  c_free(l);
  csc_spfree(A);
  c_free(q);
  csc_spfree(P);
  
  // Cost matrices
  csc_spfree(R);
  c_free(QN);
  csc_spfree(Q);

  // System matrices
  c_free(Bd);
  c_free(Ad);
  
  return 0;
}
