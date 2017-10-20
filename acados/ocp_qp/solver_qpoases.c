#include "acados/ocp_qp/ocp_qp_condensing_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/solver_qpoases.h"

#include <stdlib.h>

#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"

/* Ignore compiler warnings from qpOASES */
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-pointer-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblemB.h"
#include "qpOASES_e/QProblem.h"
#pragma clang diagnostic pop
#elif defined(__GNUC__)
    #if __GNUC__ >= 6
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
        #pragma GCC diagnostic pop
    #else
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
    #endif
#else
    #include "qpOASES_e/QProblemB.h"
    #include "qpOASES_e/QProblem.h"
#endif

#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_cond.h"

int solver_qpoases_calculate_args_size(const dense_qp_in *qp_in) {
    int size = 0;
    size += sizeof(solver_qpoases_args);
    return size;
}

char *solver_qpoases_assign_args(const dense_qp_in *qp_in,
    solver_qpoases_args **args, void *mem) {

    char *c_ptr = (char *) mem;

    *args = (solver_qpoases_args *) c_ptr;
    c_ptr += sizeof(solver_qpoases_args);

    return c_ptr;
}

static void solver_qpoases_initialize_default_args(const dense_qp_in *qp_in,
    solver_qpoases_args *args) {

    args->cputime = 1000.0;  // maximum cpu time in seconds
    args->warm_start = 0;
    args->nwsr = 1000;

}

solver_qpoases_args *solver_qpoases_create_arguments(const dense_qp_in *qp_in) {
    void *mem = malloc(solver_qpoases_calculate_args_size(qp_in));
    solver_qpoases_args *args;
    solver_qpoases_assign_args(qp_in, &args, mem);
    solver_qpoases_initialize_default_args(qp_in, args);

    return args;
}

int solver_qpoases_calculate_workspace_size(const dense_qp_in *qp_in, void *args_) {
    return 0;
}

int solver_qpoases_calculate_memory_size(const dense_qp_in *qp_in, void *args_) {
    // solver_qpoases_args *args = (solver_qpoases_args *) args_;

    int nvd = qp_in->nvd;
    // int nbd = (int)qp_in->nb;
    int ngd = qp_in->ngd;
    int ned = 0;

    // size in bytes
    int size = sizeof(solver_qpoases_memory);

    // size += 1 * nvd * nvd * sizeof(double);        // H_rmj
    size += 1 * ngd * nvd * sizeof(double);        // C_rmj
    size += 1 * ned * nvd * sizeof(double);        // A_rmj
    size += 1 * ned * sizeof(double);              // b
    size += 2 * nvd * sizeof(double);              // full_lb full_ub
    
    if (ngd > 0)  // QProblem
        size += sizeof(QProblem);
    else  // QProblemB
        size += sizeof(QProblemB);

    size = (size + 63) / 64 * 64;  // make multipl of typical cache line size
    size += 1 * 64;                // align once to typical cache line size

    return size;
}

char *solver_qpoases_assign_memory(const dense_qp_in *qp_in, void *args_,
    void **mem, void *raw_memory) {

char *c_ptr_QPdata;
// solver_qpoases_args *args = (solver_qpoases_args *) args_;
solver_qpoases_memory **qpoases_memory = (solver_qpoases_memory **) mem;

// extract problem size
int nvd = (int)qp_in->nvd;
// int nbd = (int)qp_in->nb;
int ngd = (int)qp_in->ngd;
int ned = 0;

// char pointer
char *c_ptr = (char *)raw_memory;

*qpoases_memory = (solver_qpoases_memory *) c_ptr;
c_ptr += sizeof(solver_qpoases_memory);

// align memory to typical cache line size
size_t s_ptr = (size_t)c_ptr;
s_ptr = (s_ptr + 63) / 64 * 64;
c_ptr = (char *)s_ptr;

// after alignment
c_ptr_QPdata = c_ptr;

// double stuff
// (*qpoases_memory)->H_rmj = (double *)c_ptr;
// c_ptr += nvd * nvd * sizeof(double);

(*qpoases_memory)->C_rmj = (double *)c_ptr;
c_ptr += ngd * nvd * sizeof(double);

(*qpoases_memory)->A_rmj = (double *)c_ptr;
c_ptr += ned * nvd * sizeof(double);

(*qpoases_memory)->b = (double *)c_ptr;
c_ptr += ned * sizeof(double);

(*qpoases_memory)->full_lb = (double *)c_ptr;
c_ptr += nvd * sizeof(double);

(*qpoases_memory)->full_ub = (double *)c_ptr;
c_ptr += nvd * sizeof(double);

// qpOASES (HUGE!!!)
//
if (ngd > 0) {  // QProblem
(*qpoases_memory)->QP = (void *)c_ptr;
c_ptr += sizeof(QProblem);
} else {  // QProblemB
(*qpoases_memory)->QPB = (void *)c_ptr;
c_ptr += sizeof(QProblemB);
}

// important for assigning memories for qpoases
for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
*idx = 0;

return c_ptr;
}

solver_qpoases_memory *solver_qpoases_create_memory(const dense_qp_in *qp_in,
    void *args_) {
solver_qpoases_args *args = (solver_qpoases_args *) args_;
solver_qpoases_memory *mem;
int_t memory_size = solver_qpoases_calculate_memory_size(qp_in, args);
void *raw_memory_ptr = malloc(memory_size);
char *ptr_end =
solver_qpoases_assign_memory(qp_in, args, (void **) &mem, raw_memory_ptr);
assert((char*)raw_memory_ptr + memory_size >= ptr_end); (void) ptr_end;
return mem;
}

int solver_qpoases(const dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_,
    void *memory_) {
// cast structures
solver_qpoases_args *args =
(solver_qpoases_args *)args_;
solver_qpoases_memory *memory =
(solver_qpoases_memory *)memory_;

// ocp_qp_condensing_memory *memory_cond =
// (ocp_qp_condensing_memory *)memory_cond_;

// initialize return code
int acados_status = ACADOS_SUCCESS;

int ii;
int jj;

// extract memory
// double *H_rmj = memory->H_rmj;
double *C_rmj = memory->C_rmj;
// double *A_rmj = memory->A_rmj;
// double *b = memory->b;
double *full_lb = memory->full_lb;
double *full_ub = memory->full_ub;
QProblemB *QPB = memory->QPB;
QProblem *QP = memory->QP;

// extract dense problem size and data
int nvd = (int)qp_in->nvd;
int nbd = (int)qp_in->nbd;
int ngd = (int)qp_in->ngd;
// int ned = 0;
double *H = (double *)qp_in->H;
double *g = (double *)qp_in->g;    
double *lb = (double *)qp_in->lb;
double *ub = (double *)qp_in->ub;
double *C = (double *)qp_in->C;
double *lc = (double *)qp_in->lc;
double *uc = (double *)qp_in->uc;
int *idxb = (int *)qp_in->idxb;

// convert by hand
// for (ii=0;ii<nvd;ii++){
//     for(jj=0;jj<nvd;jj++)
//         H_rmj[ii*nvd+jj] = H[jj*nvd+ii];
// }
for (ii=0;ii<ngd;ii++){
    for(jj=0;jj<nvd;jj++)
        C_rmj[ii*nvd+jj] = C[jj*ngd+ii];
}
// for (ii=0;ii<nbd;ii++){
//     for(jj=0;jj<nvd;jj++)
//         A_rmj[ii*nvd+jj] = A[jj*nvd+ii];
// }

// extract outputs
double *prim_sol = qp_out->u;
double *dual_sol = qp_out->lam;

// use qpd to convert
// struct d_dense_qp *qpd = memory_cond->qpd;
// d_cvt_dense_qp_to_rowmaj(qpd, H_rmj, g, A_rmj, b, idxb, lb, ub, C_rmj, lc, uc,
    // NULL, NULL, NULL, NULL, NULL);


// reorder bounds
for (ii = 0; ii < nvd; ii++) {
    full_lb[ii] = -QPOASES_INFTY;
    full_ub[ii] = +QPOASES_INFTY;
}
for (ii = 0; ii < nbd; ii++) {
    full_lb[idxb[ii]] = lb[ii];
    full_ub[idxb[ii]] = ub[ii];
}

// cold start the dual solution with no active constraints
int warm_start = args->warm_start;
if (!warm_start) {
    for (ii = 0; ii < 2 * nvd + 2 * ngd; ii++)
        dual_sol[ii] = 0;
}

// solve dense qp
int nwsr = args->nwsr;  // max number of working set recalculations
double cputime = args->cputime;
int return_flag = 0;
if (ngd > 0) {  // QProblem
    QProblemCON(QP, nvd, ngd, HST_POSDEF);
    QProblem_setPrintLevel(QP, PL_MEDIUM);
    QProblem_printProperties(QP);
    return_flag =
        QProblem_initW(QP, H, g, C_rmj, full_lb, full_ub, lc, uc, &nwsr, &cputime,
            NULL, dual_sol, NULL, NULL, NULL);  // NULL or 0
    //            NULL, NULL, NULL, NULL);
    //            NULL, NULL, NULL, R);  // to provide Cholesky factor
    QProblem_getPrimalSolution(QP, prim_sol);
    QProblem_getDualSolution(QP, dual_sol);
} else {  // QProblemB
    QProblemBCON(QPB, nvd, HST_POSDEF);
    QProblemB_setPrintLevel(QPB, PL_MEDIUM);
    QProblemB_printProperties(QPB);
    return_flag = QProblemB_initW(QPB, H, g, full_lb, full_ub, &nwsr, &cputime,
        NULL, dual_sol, NULL, NULL);  // NULL or 0
    QProblemB_getPrimalSolution(QPB, prim_sol);
    QProblemB_getDualSolution(QPB, dual_sol);
}

// save solution statistics to memory
memory->cputime = cputime;
memory->nwsr = nwsr;

acados_status = return_flag;
return acados_status;
}