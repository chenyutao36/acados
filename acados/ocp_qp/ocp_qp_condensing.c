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

//  #include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"

 #include <stdlib.h>
 #include <assert.h>
 
 #include "acados/utils/math.h"
 
 #include "blasfeo/include/blasfeo_target.h"
 #include "blasfeo/include/blasfeo_common.h"
 #include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
 #include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
 #include "blasfeo/include/blasfeo_d_aux.h"
 #include "blasfeo/include/blasfeo_d_blas.h"

 #include "hpipm/include/hpipm_d_ocp_qp.h"
 #include "hpipm/include/hpipm_d_ocp_qp_sol.h"
 #include "hpipm/include/hpipm_d_dense_qp.h"
 #include "hpipm/include/hpipm_d_dense_qp_sol.h"
 #include "hpipm/include/hpipm_d_cond.h"
 
 #include "acados/ocp_qp/ocp_qp_common.h" 
 #include "acados/ocp_qp/ocp_qp_condensing_common.h" 
 #include "acados/ocp_qp/ocp_qp_condensing.h"
 
 int ocp_qp_condensing_calculate_args_size(const ocp_qp_in *qp_in) {
     int N = qp_in->N;
     int size = 0;
     size += sizeof(ocp_qp_condensing_args); 

    //  size = (size+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    //  size += ALIGNMENT;

     size += (N+1)*sizeof(int);
     return size;
 }

 char *ocp_qp_condensing_assign_args(const ocp_qp_in *qp_in,
    ocp_qp_condensing_args **args, void *mem) {

    int N = qp_in->N;

    char *c_ptr = (char *) mem;

    *args = (ocp_qp_condensing_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_args);

    (*args)->scrapspace = c_ptr;
    c_ptr += (N+1)*sizeof(int);

    return c_ptr;
}

static void ocp_qp_condensing_initialize_default_args(const ocp_qp_in *qp_in,
    ocp_qp_condensing_args *args) {

    int N = qp_in->N;

    int *ns = (int *) args->scrapspace;
    int ii;
    for (ii=0; ii < N+1; ii++) ns[ii] = 0;
}



ocp_qp_condensing_args *ocp_qp_condensing_create_arguments(const ocp_qp_in *qp_in) {
    void *mem = malloc(ocp_qp_condensing_calculate_args_size(qp_in));
    ocp_qp_condensing_args *args;
    ocp_qp_condensing_assign_args(qp_in, &args, mem);
    ocp_qp_condensing_initialize_default_args(qp_in, args);

    return args;
}


int ocp_qp_condensing_calculate_workspace_size(const ocp_qp_in *qp_in, void *args_) {
    return 0;
}

int ocp_qp_condensing_calculate_memory_size(const ocp_qp_in *qp_in, void *args_) {
    ocp_qp_condensing_args *args = (ocp_qp_condensing_args *) args_;
    // extract ocp_qp_in size
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;
    int **hidxb = (int **)qp_in->idxb;

    // extract ns from args
    int_t *ns = (int_t *) args->scrapspace;

    // dummy ocp qp
    struct d_ocp_qp qp;
    qp.N = N;
    qp.nx = nx;
    qp.nu = nu;
    qp.nb = nb;
    qp.ng = ng;
    qp.idxb = hidxb;
    qp.ns = ns;

    // compute dense qp size
    int nvd = 0;
    int ned = 0;
    int nbd = 0;
    int ngd = 0;
    int_t nsd = 0;
#if 0
    // [u; x] order
    d_compute_qp_size_ocp2dense(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#else
    // [x; u] order  // XXX update with ns !!!!!
    d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#endif

    // dummy dense qp
    struct d_dense_qp qpd;
    qpd.nv = nvd;
    qpd.ne = ned;
    qpd.nb = nbd;
    qpd.ng = ngd;
    qpd.ns = nsd;

    // size in bytes
    int size = sizeof(ocp_qp_condensing_memory);

    size += 1 * sizeof(struct d_ocp_qp);                       // qp
    size += 1 * sizeof(struct d_dense_qp);                     // qpd
    size += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);  // cond_workspace

    size += 1 * sizeof(struct d_ocp_qp_sol);
    size += 1 * sizeof(struct d_dense_qp_sol);

    size += d_memsize_ocp_qp(N, nx, nu, nb, ng, ns);
    size += d_memsize_dense_qp(nvd, ned, nbd, ngd, nsd);
    size += d_memsize_cond_qp_ocp2dense(&qp, &qpd);

    size += d_memsize_ocp_qp_sol(N, nx, nu, nb, ng, ns);
    size += d_memsize_dense_qp_sol(nvd, ned, nbd, ngd, nsd);

    size += 4 * (N + 1) * sizeof(double *);  // lam_lb lam_ub lam_lg lam_ug

    size += 1 * (N + 1) * sizeof(int *);  // hidxb_rev
    

    for (int ii = 0; ii <= N; ii++) {
        size += nb[ii]*sizeof(int);  // hidxb_rev
    }  

    size = (size + 63) / 64 * 64;  // make multipl of typical cache line size
    size += 1 * 64;                // align once to typical cache line size

    return size;
}

char *ocp_qp_condensing_assign_memory(const ocp_qp_in *qp_in, void *args_,
    void **mem, void *raw_memory) {

ocp_qp_condensing_args *args = (ocp_qp_condensing_args *) args_;
ocp_qp_condensing_memory **memory = (ocp_qp_condensing_memory **) mem;

// extract problem size
int N = qp_in->N;
int *nx = (int *)qp_in->nx;
int *nu = (int *)qp_in->nu;
int *nb = (int *)qp_in->nb;
int *ng = (int *)qp_in->nc;
int **hidxb = (int **)qp_in->idxb;

// extract ns from args
int_t *ns = (int_t *) args->scrapspace;

// compute dense qp size
int nvd = 0;
int ned = 0;
int nbd = 0;
int ngd = 0;
int_t nsd = 0;
#if 0
// [u; x] order
d_compute_qp_size_ocp2dense(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#else
// [x; u] order  // XXX update with ns !!!!!
d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#endif


// char pointer
char *c_ptr = (char *)raw_memory;

*memory = (ocp_qp_condensing_memory *) c_ptr;
c_ptr += sizeof(ocp_qp_condensing_memory);

(*memory)->qp = (struct d_ocp_qp *)c_ptr;
c_ptr += 1 * sizeof(struct d_ocp_qp);

(*memory)->qpd = (struct d_dense_qp *)c_ptr;
c_ptr += 1 * sizeof(struct d_dense_qp);

(*memory)->cond_workspace =
(struct d_cond_qp_ocp2dense_workspace *)c_ptr;
c_ptr += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);

(*memory)->qp_sol = (struct d_ocp_qp_sol *)c_ptr;
c_ptr += 1 * sizeof(struct d_ocp_qp_sol);

(*memory)->qpd_sol = (struct d_dense_qp_sol *)c_ptr;
c_ptr += 1 * sizeof(struct d_dense_qp_sol);

(*memory)->hlam_lb = (double **)c_ptr;
c_ptr += (N + 1) * sizeof(double *);
//
(*memory)->hlam_ub = (double **)c_ptr;
c_ptr += (N + 1) * sizeof(double *);
//
(*memory)->hlam_lg = (double **)c_ptr;
c_ptr += (N + 1) * sizeof(double *);
//
(*memory)->hlam_ug = (double **)c_ptr;
c_ptr += (N + 1) * sizeof(double *);

(*memory)->hidxb_rev = (int **)c_ptr;
c_ptr += (N + 1) * sizeof(int *);

struct d_ocp_qp *qp = (*memory)->qp;

struct d_dense_qp *qpd = (*memory)->qpd;

struct d_cond_qp_ocp2dense_workspace *cond_workspace =
(*memory)->cond_workspace;

struct d_ocp_qp_sol *qp_sol = (*memory)->qp_sol;

struct d_dense_qp_sol *qpd_sol = (*memory)->qpd_sol;

// align memory to typical cache line size
size_t s_ptr = (size_t)c_ptr;
s_ptr = (s_ptr + 63) / 64 * 64;
c_ptr = (char *)s_ptr;

// ocp qp structure
d_create_ocp_qp(N, nx, nu, nb, ng, ns, qp, c_ptr);
c_ptr += qp->memsize;
// dense qp structure
d_create_dense_qp(nvd, ned, nbd, ngd, nsd, qpd, c_ptr);
c_ptr += qpd->memsize;
// cond workspace structure
d_create_cond_qp_ocp2dense(qp, qpd, cond_workspace, c_ptr);
c_ptr += cond_workspace->memsize;

d_create_ocp_qp_sol(N, nx, nu, nb, ng, ns, qp_sol, c_ptr);
c_ptr += qp_sol->memsize;

d_create_dense_qp_sol(nvd, ned, nbd, ngd, nsd, qpd_sol, c_ptr);
c_ptr += qpd_sol->memsize;


for (int ii = 0; ii <= N; ii++) {
(*memory)->hidxb_rev[ii] = (int *) c_ptr;
c_ptr += nb[ii]*sizeof(int);
}

return c_ptr;
}

ocp_qp_condensing_memory *ocp_qp_condensing_create_memory(const ocp_qp_in *qp_in,
    void *args_) {
ocp_qp_condensing_args *args = (ocp_qp_condensing_args *) args_;
ocp_qp_condensing_memory *mem;
int_t memory_size = ocp_qp_condensing_calculate_memory_size(qp_in, args);
void *raw_memory_ptr = malloc(memory_size);
char *ptr_end =
ocp_qp_condensing_assign_memory(qp_in, args, (void **) &mem, raw_memory_ptr);
assert((char*)raw_memory_ptr + memory_size >= ptr_end); (void) ptr_end;
return mem;
}


int ocp_qp_condensing(const ocp_qp_in *qp_in, dense_qp_in *qp_out
    ,void *memory_) {
// cast structures

ocp_qp_condensing_memory *memory =
(ocp_qp_condensing_memory *)memory_; 

// loop index
int ii, jj;

// extract memory
struct d_ocp_qp *qp = memory->qp;
struct d_dense_qp *qpd = memory->qpd;
struct d_cond_qp_ocp2dense_workspace *cond_workspace =
memory->cond_workspace;

int **hidxb_rev = (int **) memory->hidxb_rev;

// extract ocp problem size
int N = qp_in->N;
int *nx = (int *) qp_in->nx;
int *nu = (int *) qp_in->nu;
int *nb = (int *) qp_in->nb;
// int *ng = (int *) qp_in->nc;

// extract input data
double **hA = (double **)qp_in->A;
double **hB = (double **)qp_in->B;
double **hb = (double **)qp_in->b;
double **hQ = (double **)qp_in->Q;
double **hS = (double **)qp_in->S;
double **hR = (double **)qp_in->R;
double **hq = (double **)qp_in->q;
double **hr = (double **)qp_in->r;
double **hd_lb = (double **)qp_in->lb;
double **hd_ub = (double **)qp_in->ub;
double **hC = (double **)qp_in->Cx;
double **hD = (double **)qp_in->Cu;
double **hd_lg = (double **)qp_in->lc;
double **hd_ug = (double **)qp_in->uc;
int **hidxb = (int **)qp_in->idxb;

// extract output properties
int nvd = qp_out->nvd;
double *H = qp_out->H;
double *g = qp_out->g;
double *A = qp_out->A;
double *b = qp_out->b;
double *d_lb0 = qp_out->lb;
double *d_ub0 = qp_out->ub;
double *d_lg = qp_out->lc;
double *d_ug = qp_out->uc;
double *C = qp_out->C;
int *idxb = qp_out->idxb;

// compute bounds indeces in order [u; x]
for (ii = 0; ii <= N; ii++) {
    for (jj = 0; jj < nb[ii]; jj++) {
        if (hidxb[ii][jj] < nx[ii]) {  // state constraint
            hidxb_rev[ii][jj] = hidxb[ii][jj]+nu[ii];
        } else  {  // input constraint
            hidxb_rev[ii][jj] = hidxb[ii][jj]-nx[ii];
        }
    }
}

// ocp qp structure
d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb_rev, hd_lb, hd_ub, hC, hD,
    hd_lg, hd_ug, NULL, NULL, NULL, NULL, NULL, qp);

// dense qp structure
d_cond_qp_ocp2dense(qp, qpd, cond_workspace);
// d_print_tran_strvec(nvd, qpd->g, 0);

// fill in the upper triangular of H in dense_qp
dtrtr_l_libstr(nvd, qpd->Hg, 0, 0, qpd->Hg, 0, 0);
// d_print_tran_strvec(nvd, qpd->g, 0);

// convert the qpd struct into dense_qp_in class
d_cvt_dense_qp_to_colmaj(qpd, H, g, A, b, idxb, d_lb0, d_ub0, C, d_lg, d_ug,
    NULL, NULL, NULL, NULL, NULL);
// printf("\nin condensing\n");
// d_print_tran_strvec(nvd, qpd->g, 0);
// d_print_mat(1, nvd, g, 1);

// return
return 0;
//
}

void ocp_qp_condensing_initialize(const ocp_qp_in *qp_in, void *args_, void **mem,
    void **work) {
ocp_qp_condensing_args *args = (ocp_qp_condensing_args *)args_;

*mem = ocp_qp_condensing_create_memory(qp_in, args);
*work = NULL;
}

void ocp_qp_condensing_destroy(void *mem_, void *work) {
    (void) mem_;
    (void) work;
}