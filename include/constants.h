#ifndef CONSTANTS_H
#define CONSTANTS_H


/*******************
 * OSQP Versioning *
 *******************/
#define OSQP_VERSION ("0.1") /* string literals automatically null-terminated */



/******************
 * Solver Status  *
 ******************/
 // TODO: Add other statuses
#define OSQP_SOLVED (1)
// #define OSQP_SOLVED_INACCURATE (2)
#define OSQP_INFEASIBLE (-2) /* primal infeasible, dual unbounded   */
#define OSQP_UNSOLVED (-10)  /* Unsolved. Only setup function has been called */



/**********************************
 * Solver Parameters and Settings *
 **********************************/
/* TODO: complete parameters, these are just a couple of examples */

#define NORMALIZE (0)
#define RHO (1.6)
#define MAX_ITER (2500)
#define EPS_ABS (1E-5)
#define EPS_REL (1E-5)
#define ALPHA (1.6)
#define VERBOSE (1)
#define WARM_START (0)
#define DELTA (1E-7)


/* Printing */
#define PRINT_INTERVAL 10
















#endif