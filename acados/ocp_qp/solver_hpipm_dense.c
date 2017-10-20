#include "acados/ocp_qp/ocp_qp_condensing_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/solver_hpipm_dense.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"

#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

int solver_hpipm_calculate_args_size(const dense_qp_in *qp_in) {
    int size = 0;
    size += sizeof(ocp_qp_condensing_hpipm_args);
    return size;
}

char *solver_hpipm_assign_args(const dense_qp_in *qp_in,
    solver_hpipm_args **args, void *mem) {


    char *c_ptr = (char *) mem;

    *args = (ocp_qp_condensing_hpipm_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_hpipm_args);

    return c_ptr;
}

static void solver_hpipm_initialize_default_args(const dense_qp_in *qp_in,
    solver_hpipm_args *args) {
    
    args->res_g_max = 1e-6;
    args->res_b_max = 1e-8;
    args->res_d_max = 1e-8;
    args->res_m_max = 1e-8;
    args->iter_max = 50;
    args->alpha_min = 1e-8;
    args->mu0 = 1;
}

solver_hpipm_args *solver_hpipm_create_arguments(const dense_qp_in *qp_in) {
    void *mem = malloc(solver_hpipm_calculate_args_size(qp_in));
    solver_hpipm_args *args;
    solver_hpipm_assign_args(qp_in, &args, mem);
    solver_hpipm_initialize_default_args(qp_in, args);

    return args;
}

int_t solver_hpipm_calculate_workspace_size(const dense_qp_in *qp_in, void *args_) {
    return 0;
}

int_t solver_hpipm_calculate_memory_size(const dense_qp_in *qp_in, void *args_) {
    solver_hpipm_args *args = (solver_hpipm_args *) args_;
    // extract ocp qp in size
    int nvd = qp_in->nvd;
    int nbd = qp_in->nb;
    int ngd = qp_in->ngd;
    int ned = 0;
    int nsd = 0;

    // dummy dense qp
    struct d_dense_qp qpd;
    qpd.nv = nvd;
    qpd.ne = ned;
    qpd.nb = nbd;
    qpd.ng = ngd;
    qpd.ns = nsd;

    // dummy ipm arg
    struct d_dense_qp_ipm_arg ipm_arg;
    ipm_arg.stat_max = args->iter_max;

    // size in bytes
    int_t size = sizeof(solver_hpipm_memory);

    size += 1 * sizeof(struct d_dense_qp);                     // qpd
    size += 1 * sizeof(struct d_dense_qp_sol);                 // qpd_sol
    size += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);  // cond_workspace
    size += 1 * sizeof(struct d_dense_qp_ipm_arg);        // ipm_arg
    size += 1 * sizeof(struct d_dense_qp_ipm_workspace);  // ipm_workspace

    size += d_memsize_dense_qp(nvd, ned, nbd, ngd, nsd);
    size += d_memsize_dense_qp_sol(nvd, ned, nbd, ngd, nsd);
    size += d_memsize_cond_qp_ocp2dense(&qp, &qpd);
    size += d_memsize_dense_qp_ipm_arg(&qpd);
    size += d_memsize_dense_qp_ipm(&qpd, &ipm_arg);

    size = (size + 63) / 64 * 64;  // make multipl of typical cache line size
    size += 1 * 64;                // align once to typical cache line size

    return size;
}

char *solver_hpipm_assign_memory(const dense_qp_in *qp_in, void *args_, void **mem_,
    void *raw_memory) {

solver_hpipm_args *args = (solver_hpipm_args *) args_;
solver_hpipm_memory **hpipm_memory = (solver_hpipm_memory **) mem_;
// extract ocp qp in size
int nvd = qp_in->nvd;
int nbd = qp_in->nb;
int ngd = qp_in->ngd;
int ned = 0;
int nsd = 0;

// char pointer
char *c_ptr = (char *)raw_memory;

*hpipm_memory = (ocp_qp_condensing_hpipm_memory *) c_ptr;
c_ptr += sizeof(ocp_qp_condensing_hpipm_memory);

//
(*hpipm_memory)->qpd = (struct d_dense_qp *)c_ptr;
c_ptr += 1 * sizeof(struct d_dense_qp);
//
(*hpipm_memory)->qpd_sol = (struct d_dense_qp_sol *)c_ptr;
c_ptr += 1 * sizeof(struct d_dense_qp_sol);
//
(*hpipm_memory)->cond_workspace =
(struct d_cond_qp_ocp2dense_workspace *)c_ptr;
c_ptr += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);
//
(*hpipm_memory)->ipm_arg = (struct d_dense_qp_ipm_arg *)c_ptr;
c_ptr += 1 * sizeof(struct d_dense_qp_ipm_arg);
//
(*hpipm_memory)->ipm_workspace = (struct d_dense_qp_ipm_workspace *)c_ptr;
c_ptr += 1 * sizeof(struct d_dense_qp_ipm_workspace);

struct d_dense_qp *qpd = (*hpipm_memory)->qpd;
//
struct d_dense_qp_sol *qpd_sol = (*hpipm_memory)->qpd_sol;
//
struct d_cond_qp_ocp2dense_workspace *cond_workspace =
(*hpipm_memory)->cond_workspace;
//
struct d_dense_qp_ipm_arg *ipm_arg = (*hpipm_memory)->ipm_arg;
//
struct d_dense_qp_ipm_workspace *ipm_workspace =
(*hpipm_memory)->ipm_workspace;

// align memory to typical cache line size
size_t s_ptr = (size_t)c_ptr;
s_ptr = (s_ptr + 63) / 64 * 64;
c_ptr = (char *)s_ptr;


// dense qp structure
d_create_dense_qp(nvd, ned, nbd, ngd, nsd, qpd, c_ptr);
c_ptr += qpd->memsize;
// dense qp sol structure
d_create_dense_qp_sol(nvd, ned, nbd, ngd, nsd, qpd_sol, c_ptr);
c_ptr += qpd_sol->memsize;
// cond workspace structure
d_create_cond_qp_ocp2dense(qp, qpd, cond_workspace, c_ptr);
c_ptr += cond_workspace->memsize;
// ipm arg structure
d_create_dense_qp_ipm_arg(qpd, ipm_arg, c_ptr);
c_ptr += ipm_arg->memsize;
d_set_default_dense_qp_ipm_arg(ipm_arg);
ipm_arg->iter_max = args->iter_max;
ipm_arg->stat_max = args->iter_max;
ipm_arg->alpha_min = args->alpha_min;
ipm_arg->res_g_max = args->res_g_max;
ipm_arg->res_b_max = args->res_b_max;
ipm_arg->res_d_max = args->res_d_max;
ipm_arg->res_m_max = args->res_m_max;
ipm_arg->mu0 = args->mu0;
// ipm workspace structure
d_create_dense_qp_ipm(qpd, ipm_arg, ipm_workspace, c_ptr);
c_ptr += ipm_workspace->memsize;

return c_ptr;
}

solver_hpipm_memory *solver_hpipm_create_memory(const dense_qp_in *qp_in,
    void *args_) {

solver_hpipm_args *args = (solver_hpipm_args *) args_;

solver_hpipm_memory *mem;
int_t memory_size = solver_hpipm_calculate_memory_size(qp_in, args);
void *raw_memory = malloc(memory_size);
char *ptr_end =
solver_hpipm_assign_memory(qp_in, args, (void **) &mem, raw_memory);
assert((char*) raw_memory + memory_size >= ptr_end); (void) ptr_end;

return mem;
}

int_t solver_hpipm(const dense_qp_in *qp_in, dense_qp_out *qp_out,
    void *args_,
    void *memory_,
    void *workspace_) {

    solver_hpipm_args *args = (solver_hpipm_args *) args_;
    solver_hpipm_memory *memory = (solver_hpipm_memory *) memory_;

// initialize return code
int_t acados_status = ACADOS_SUCCESS;

// loop index
int_t ii, jj;

struct d_dense_qp *qpd = memory->qpd;
struct d_dense_qp_sol *qpd_sol = memory->qpd_sol;
struct d_cond_qp_ocp2dense_workspace *cond_workspace =
memory->cond_workspace;
struct d_dense_qp_ipm_arg *ipm_arg = memory->ipm_arg;
struct d_dense_qp_ipm_workspace *ipm_workspace = memory->ipm_workspace;

// extract problem size
int nvd = qp_in->nvd;
int nbd = qp_in->nb;
int ngd = qp_in->ngd;
int ned = 0;
int nsd = 0;

double *H = (double *)qp_in->H;
double *g = (double *)qp_in->g; 
double *A = (double *)qp_in->A;
double *b = (double *)qp_in->b;      
double *lb = (double *)qp_in->lb;
double *ub = (double *)qp_in->ub;
double *C = (double *)qp_in->C;
double *lc = (double *)qp_in->lc;
double *uc = (double *)qp_in->uc;
int *idxb = (int *)qp_in->idxb;

double *prim_sol = qp_out->u;
double *dual_sol = qp_out->lam;

// convert col_maj to qpd struct
d_cvt_colmaj_to_dense_qp(H, g, A, b, idxb, lb, ub, C, lc, uc, 
    NULL, NULL, NULL, NULL, NULL, qpd);

// solve ipm
d_solve_dense_qp_ipm(qpd, qpd_sol, ipm_arg, ipm_workspace);


// convert dense_qp_sol to col_maj
d_cvt_dense_qp_sol_to_colmaj(qpd, qpd_sol, v, ls, us, pi, lam_lb, lam_ub, lam_lg, lam_ug, lam_ls, lam_us);

// extract iteration number
memory->iter = ipm_workspace->iter;

// compute infinity norm of residuals
real_t *inf_norm_res = memory->inf_norm_res;
real_t res_tmp;
//
real_t *res_g;
inf_norm_res[0] = 0;
res_g = ipm_workspace->res_g->pa;
for (jj = 0; jj < qpd->nv; jj++) {
res_tmp = fabs(res_g[jj]);
if (res_tmp > inf_norm_res[0]) {
inf_norm_res[0] = res_tmp;
}
}
real_t *res_b;
inf_norm_res[1] = 0;
res_b = ipm_workspace->res_b->pa;
for (jj = 0; jj < qpd->ne; jj++) {
res_tmp = fabs(res_b[jj]);
if (res_tmp > inf_norm_res[1]) {
inf_norm_res[1] = res_tmp;
}
}
real_t *res_d;
inf_norm_res[2] = 0;
res_d = ipm_workspace->res_d->pa;
for (jj = 0; jj < 2 * qpd->nb + 2 * qpd->ng; jj++) {
res_tmp = fabs(res_d[jj]);
if (res_tmp > inf_norm_res[2]) {
inf_norm_res[2] = res_tmp;
}
}
real_t *res_m;
inf_norm_res[3] = 0;
res_m = ipm_workspace->res_m->pa;
for (jj = 0; jj < 2 * qpd->nb + 2 * qpd->ng; jj++) {
res_tmp = fabs(res_m[jj]);
if (res_tmp > inf_norm_res[3]) {
inf_norm_res[3] = res_tmp;
}
}
inf_norm_res[4] = ipm_workspace->res_mu;

// max number of iterations
if (ipm_workspace->iter == args->iter_max) acados_status = ACADOS_MAXITER;
// minimum step length
if (ipm_workspace->stat[3 + (ipm_workspace->iter - 1) * 5] <
args->alpha_min)
acados_status = ACADOS_MINSTEP;

// return
return acados_status;
//
}