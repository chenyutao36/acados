#ifndef ACADOS_OCP_QP_SOLVER_QPOASES_H_
#define ACADOS_OCP_QP_SOLVER_QPOASES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_condensing_common.h"
#include "acados/utils/types.h"
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"

typedef struct solver_qpoases_args_ {
    double cputime;  // maximum cpu time in seconds
    int nwsr;        // maximum number of working set recalculations
    int warm_start;  // warm start with dual_sol in memory
} solver_qpoases_args;

typedef struct solver_qpoases_memory_ {
    // double *H_rmj;
    double *C_rmj;
    double *A_rmj;
    double *b;
    double *full_lb;
    double *full_ub;
    void *QPB;  // XXX cast to QProblemB to use !!!
    void *QP;   // XXX cast to QProblem to use !!!
    double inf_norm_res[5];
    double cputime;  // required cpu time
    int nwsr;        // performed number of working set recalculations
} solver_qpoases_memory;

typedef struct {
    int_t (*fun)(const dense_qp_in *qp_in, dense_qp_out *qp_out, void *args, void *mem, void *work);
    void (*initialize)(const dense_qp_in *qp_in, void *args, void **mem, void **work);
    void (*destroy)(void *mem, void *work);
    dense_qp_in *qp_in;
    dense_qp_out *qp_out;
    void *args;
    void *mem;
    void *work;
} dense_qp_solver;

solver_qpoases_args *solver_qpoases_create_arguments(const dense_qp_in *qp_in);

int_t solver_qpoases_calculate_memory_size(const dense_qp_in *qp_in, void *args_);

char *solver_qpoases_assign_memory(const dense_qp_in *qp_in, void *args_,
                                void **qpoases_memory, void *raw_memory);

solver_qpoases_memory *solver_qpoases_create_memory(const dense_qp_in *qp_in,
                                                                          void *args_);
int_t solver_qpoases_calculate_workspace_size(const dense_qp_in *qp_in, void *args_);

int_t solver_qpoases(const dense_qp_in *input, dense_qp_out *output, void *args_,
                    void *memory_);

void solver_qpoases_initialize(const dense_qp_in *qp_in, void *args_, void **mem,
                                          void **work);

void solver_qpoases_destroy(void *mem, void *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_SOLVER_QPOASES_H_