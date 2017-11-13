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

 #ifndef ACADOS_OCP_QP_OCP_QP_CONDENSING_COMMON_H_
 #define ACADOS_OCP_QP_OCP_QP_CONDENSING_COMMON_H_
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 #include "acados/utils/types.h"

 
 // define a dense_qp_in struct
 typedef struct{    
     int_t nvd;
     int_t nbd;
     int_t ngd;
     real_t *H;
     real_t *g;
     real_t *A;
     real_t *b;
     real_t *lb;
     real_t *ub;    
     real_t *C;
     real_t *lc;
     real_t *uc;
     int_t *idxb;
 }dense_qp_in;

typedef struct {
    int nvd;
    int nbd;
    int ngd;
    real_t *u;
    real_t *lam;
} dense_qp_out;

 
int_t dense_qp_in_calculate_size(int_t N, int_t *nx, int_t *nu, 
    int_t *nb, int_t *nc, int_t **idxb);
 
char *assign_dense_qp_in(int_t N, int_t *nx, int_t *nu, int_t *nb,
    int_t *nc, int_t **idxb, dense_qp_in **qp_in, void *ptr);
 
dense_qp_in *create_dense_qp_in(int_t N, int_t *nx, int_t *nu, 
    int_t *nb, int_t *nc, int_t **idxb);


int_t dense_qp_out_calculate_size(int_t nvd, int_t nbd, int_t ngd);
     
char *assign_dense_qp_out(int_t nvd, int_t nbd, int_t ngd, dense_qp_out **qp_out, void *ptr);
     
dense_qp_out *create_dense_qp_out(int nvd, int_t nbd, int ngd);
 
 #ifdef __cplusplus
 } /* extern "C" */
 #endif
 
 #endif  // ACADOS_OCP_QP_OCP_QP_COMMON_H_
 