#ifndef ACADOS_OCP_QP_SOLVER_HPIPM_DENSE_H_
#define ACADOS_OCP_QP_SOLVER_HPIPM_DENSE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef struct solver_hpipm_args_ {
    real_t alpha_min;
    real_t res_g_max;
    real_t res_b_max;
    real_t res_d_max;
    real_t res_m_max;
    real_t mu0;
    int_t iter_max;
} solver_hpipm_args;

typedef struct solver_hpipm_memory_ {
    struct d_dense_qp *qpd;
    struct d_dense_qp_sol *qpd_sol;
    struct d_cond_qp_ocp2dense_workspace *cond_workspace;
    struct d_dense_qp_ipm_arg *ipm_arg;
    struct d_dense_qp_ipm_workspace *ipm_workspace;
    real_t inf_norm_res[5];
    int_t iter;
} solver_hpipm_memory;

solver_hpipm_args *solver_hpipm_create_arguments(const dense_qp_in *qp_in);

int_t solver_hpipm_calculate_memory_size(const dense_qp_in *qp_in,
     void *args_);

char *solver_hpipm_assign_memory(const dense_qp_in *qp_in, void *args_, 
    void **mem_, void *raw_memory);
//
solver_hpipm_memory *solver_hpipm_create_memory(const dense_qp_in *qp_in,
    void *args_);

int_t solver_hpipm_calculate_workspace_size(const dense_qp_in *qp_in, 
    void *args_);

int_t solver_hpipm(const dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_,
    void *memory_, void *workspace_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_SOLVER_HPIPM_DENSE_H_