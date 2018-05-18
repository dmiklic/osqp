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
csc* mpc_to_osqp_P(csc* Q, csc* QN, csc* R, c_int N);

/* 
 * Creates the OSQP linear cost vector q from MPC cost matrices
 * and setpoint value.
 * 
 * Assumes constant setpoint throughout the prediction horizon.
 */
c_float* mpc_to_osqp_q_const_xr(csc* Q, csc* QN, c_float* xr, c_int nu, c_int N);

/* 
 * Creates the OSQP linear cost vector q from MPC cost matrices
 * and setpoint value.
 * 
 * xr must be an array of Q->m \times (N+1) elements, specifying setponts for each
 * step of the prediction horizon (0 through N).
 */
c_float* mpc_to_osqp_q(csc* Q, csc* QN, c_float* xr, c_int nu, c_int N);

/* 
 * Creates the OSQP constraint matrix A from system dynamics matrices Ad and Bd
 */
csc* mpc_to_osqp_A(csc* Ad, csc* Bd, c_int N);

/* 
 * Creates the OSQP bound vector from state and input bounds.
 * Provide -x0, x_min, xN_min and u_min for lower bound and
 * x0, x_max, xN_max and u_max for upper bound.
 */
c_float* mpc_to_osqp_bound(c_float* minus_x0, c_float* x_bound, c_float* xN_bound,
                           c_float* u_bound, c_int nx, c_int nu, c_int N);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif
