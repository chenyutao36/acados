/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

 #ifndef ACADOS_OCP_QP_OCP_QP_CONDENSING_H_
 #define ACADOS_OCP_QP_OCP_QP_CONDENSING_H_

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 #include "blasfeo/include/blasfeo_target.h"
 #include "blasfeo/include/blasfeo_common.h"
 #include "acados/ocp_qp/ocp_qp_common.h"
 #include "acados/ocp_qp/ocp_qp_condensing_common.h"
 #include "acados/utils/types.h"
 

typedef struct ocp_qp_condensing_args_ {
    void *scrapspace;
} ocp_qp_condensing_args;

typedef struct ocp_qp_condensing_memory_ {
    struct d_ocp_qp *qp;
    struct d_dense_qp *qpd;
    struct d_cond_qp_ocp2dense_workspace *cond_workspace;
    struct d_ocp_qp_sol *qp_sol;
    struct d_dense_qp_sol *qpd_sol;
    double **hlam_lb;
    double **hlam_ub;
    double **hlam_lg;
    double **hlam_ug;
    int **hidxb_rev;
} ocp_qp_condensing_memory;

ocp_qp_condensing_args *ocp_qp_condensing_create_arguments(const ocp_qp_in *qp_in);

int_t ocp_qp_condensing_calculate_memory_size(const ocp_qp_in *qp_in, void *args_);

char *ocp_qp_condensing_assign_memory(const ocp_qp_in *qp_in, void *args_, void **memory,
     void *raw_memory);

ocp_qp_condensing_memory *ocp_qp_condensing_create_memory(const ocp_qp_in *qp_in,
                                void *args_);

int_t ocp_qp_condensing_calculate_workspace_size(const ocp_qp_in *qp_in, void *args_);

int_t ocp_qp_condensing(const ocp_qp_in *input, dense_qp_in *output, 
void *memory_);

void ocp_qp_condensing_initialize(const ocp_qp_in *qp_in, void *args_, void **mem,
void **work);

void ocp_qp_condensing_destroy(void *mem, void *work);

int ocp_qp_expand(const ocp_qp_in *qp_in, const dense_qp_out *sol, ocp_qp_out *qp_out, void *args_,
    void *memory_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_