#include <stdlib.h>
#include <assert.h>

#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"

#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_cond.h"

#include "acados/ocp_qp/ocp_qp_condensing_common.h"


// calculate the size of dense_qp_in
// int_t dense_qp_in_calculate_size(const int_t N, const int_t *nx, const int_t *nu, 
//     const int_t *nb, const int_t *nc, const int_t **idxb) {
int_t dense_qp_in_calculate_size(int_t N, int_t *nx, int_t *nu, 
    int_t *nb, int_t *nc, int_t **idxb) {

int_t bytes = sizeof(dense_qp_in);

int_t nvd = 0;
int_t ned = 0;
int_t nbd = 0;
int_t ngd = 0;
// int_t ns_d = 0;

d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, idxb, nc, &nvd, &ned, &nbd, &ngd);

bytes += nvd * nvd * sizeof(double); //H
bytes += nvd * sizeof(double); //g
bytes += ned * nvd * sizeof(double); //A
bytes += ned * sizeof(double); //b
bytes += 2* nbd * sizeof(double); //lb, ub
bytes += ngd * nvd * sizeof(double); //C
bytes += 2* ngd * sizeof(double); //lc, uc
bytes += 1* nbd * sizeof(int); //idxb

bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
bytes += ALIGNMENT;

return bytes;
}


// assign pointer to dense_qp_in
char *assign_dense_qp_in(int_t N, int_t *nx, int_t *nu, int_t *nb,
    int_t *nc, int_t **idxb, dense_qp_in **qp_in, void *ptr) {

// pointer to initialize QP data to zero
// char *c_ptr_QPdata;

// char pointer
char *c_ptr = (char *) ptr;

*qp_in = (dense_qp_in *) c_ptr;
c_ptr += sizeof(dense_qp_in);

int_t nvd = 0;
int_t ned = 0;
int_t nbd = 0;
int_t ngd = 0;
// int_t ns_d = 0;
d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, idxb, nc, &nvd, &ned, &nbd, &ngd);

// copy dimensions to workspace

(*qp_in)->nvd = nvd;

(*qp_in)->nbd = nbd;

(*qp_in)->ngd = ngd;

// align data
size_t l_ptr = (size_t) c_ptr;
l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
c_ptr = (char *) l_ptr;

// assign pointers to doubles
// c_ptr_QPdata = c_ptr;

// assign pointers
(*qp_in)->H = (real_t *) c_ptr;
c_ptr += nvd*nvd*sizeof(real_t);

(*qp_in)->g = (real_t *) c_ptr;
c_ptr += nvd * sizeof(real_t);

(*qp_in)->A = (real_t *) c_ptr;
c_ptr += ned*nvd*sizeof(real_t);

(*qp_in)->b = (real_t *) c_ptr;
c_ptr += ned*sizeof(real_t);

(*qp_in)->lb = (real_t *) c_ptr;
c_ptr += nbd*sizeof(real_t);

(*qp_in)->ub = (real_t *) c_ptr;
c_ptr += nbd*sizeof(real_t);

(*qp_in)->C = (real_t *) c_ptr;
c_ptr += ngd*nvd*sizeof(real_t);

(*qp_in)->lc = (real_t *) c_ptr;
c_ptr += ngd*sizeof(real_t);

(*qp_in)->uc = (real_t *) c_ptr;
c_ptr += ngd*sizeof(real_t);

(*qp_in)->idxb = (int_t *) c_ptr;
c_ptr += nbd*sizeof(int_t);

// set QP data to zero (mainly for valgrind)
// for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
// *idx = 0;

return c_ptr;
}

dense_qp_in *create_dense_qp_in(int_t N, int_t *nx, int_t *nu, 
    int_t *nb, int_t *nc, int_t **idxb) {

dense_qp_in *qp_in;

int_t bytes = dense_qp_in_calculate_size(N, nx, nu, nb, nc, idxb);

void *ptr = malloc(bytes);

char *ptr_end = assign_dense_qp_in(N, nx, nu, nb, nc, idxb, &qp_in, ptr);
assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

return qp_in;
}


int_t dense_qp_out_calculate_size(int_t nvd, int_t nbd, int_t ngd){
        int_t bytes = sizeof(dense_qp_out);

        bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
        bytes += ALIGNMENT;
        
        bytes += 1 * nvd * sizeof(double);              // u
        bytes += (2 * nvd + 2 * ngd) * sizeof(double);  // lam      
        
        return bytes;
    }
     
char *assign_dense_qp_out(int_t nvd, int_t nbd, int_t ngd, dense_qp_out **qp_out, void *ptr){
// pointer to initialize QP data to zero
char *c_ptr_QPdata;

// char pointer
char *c_ptr = (char *) ptr;


*qp_out = (dense_qp_out *) c_ptr;
c_ptr += sizeof(dense_qp_out);

(*qp_out)->nvd = nvd;
(*qp_out)->nbd = nbd;
(*qp_out)->ngd = ngd;

// align data
size_t l_ptr = (size_t) c_ptr;
l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
c_ptr = (char *) l_ptr;

c_ptr_QPdata = c_ptr;

// assign pointers to doubles
(*qp_out)->u = (double *)c_ptr;
c_ptr += nvd * sizeof(double);

(*qp_out)->lam = (double *)c_ptr;
c_ptr += (2 * nvd + 2 * ngd) * sizeof(double);

// set QP data to zero (mainly for valgrind)
for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
*idx = 0;

return c_ptr;
}
     
dense_qp_out *create_dense_qp_out(int nvd, int_t nbd, int ngd){
    dense_qp_out *qp_out;
        
    int_t bytes = dense_qp_out_calculate_size(nvd, nbd, ngd);
        
    void *ptr = malloc(bytes);
        
    char *ptr_end = assign_dense_qp_out(nvd, nbd, ngd, &qp_out, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;
        
    return qp_out;
}

