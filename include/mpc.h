/* Helper functions for casting MPC problems to QP */

#ifndef MPC_H
#define MPC_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "cs.h"

/* 
 * Creates the OSQP quadratic cost matrix Q from MPC cost matrices.
 */
csc* mpc_to_osqp_P(const csc* Q,
                   const csc* QN,
                   const csc* R,
                   c_int N);

/* 
 * Creates the OSQP linear cost vector q from MPC cost matrices
 * and setpoint value.
 * 
 * Assumes constant setpoint throughout the prediction horizon.
 */
c_float* mpc_to_osqp_q_const_xr(const csc* Q,
                                const csc* QN,
                                const c_float* xr,
                                c_int nu, c_int N);

/* 
 * Creates the OSQP linear cost vector q from MPC cost matrices
 * and setpoint value.
 * 
 * xr must be an array of Q->m \times (N+1) elements, specifying setponts for each
 * step of the prediction horizon (0 through N).
 */
c_float* mpc_to_osqp_q(csc* Q,
                       csc* QN,
                       const c_float* xr,
                       c_int nu, c_int N);

/* 
 * Creates the OSQP constraint matrix A from system dynamics matrices Ad and Bd
 */
csc* mpc_to_osqp_A(const csc* Ad,
                   const csc* Bd,
                   c_int N);

/*
 * Creates the OSQP constraint matrix A for a LPV-A-MPC problem, i.e.,
 * an LPV-MPC problem where the state transition matrix Ad changes at
 * every prediction step, while the input matrix Bd remains constant
 * over the prediction horizon.
 *
 * The array Ad[] must contain exactly N csc* objects.
 * It is assumed that Ad[k]->nzmax is the same for all k \in [0,N)
 *
 * Uineq is an array of matrices specifying additional inequality
 * constraints that must be satisfied by system inputs at each
 * prediction step.
 */
csc* lpv_a_mpc_to_osqp_A(csc *Ad[], const csc* Bd, csc *Uineq[], c_int N);

/* 
 * Creates the OSQP bound vector from state and input bounds.
 * Provide -x0, x_min, xN_min and u_min for lower bound and
 * x0, x_max, xN_max and u_max for upper bound.
 */
c_float* mpc_to_osqp_bound(const c_float* minus_x0,
                           const c_float* x_bound,
                           const c_float* xN_bound,
                           const c_float* u_bound,
                           c_int nx, c_int nu, c_int N);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif
