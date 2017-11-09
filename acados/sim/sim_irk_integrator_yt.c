// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_irk_integrator_yt.h"

#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"


int irk_calculate_memory_size(sim_RK_opts *opts, sim_in *in)
{

    int nx = in->nx;
    int nu = in->nu; 
    // int NF = in->NF;
    
    int ns = opts->num_stages; // number of stages

    int size = sizeof(sim_irk_memory);

    size += 1*sizeof(struct d_strmat); // JG
    size += 2*sizeof(struct d_strvec); // rG, K

    size += 1*d_size_strmat(nx*ns, nx*ns); // JG
    size += 2*d_size_strmat(nx*ns, 1); // rG K

    size += 2*nx*sizeof(double); // xt0, xt1
    size += nx*sizeof(double); //kt
    size += nx*sizeof(double); //  rGt
    size += 2*nx*nx*sizeof(double); // jac_out
    size += nx*nx*sizeof(double); // Jt
    size += (2*nx+nu)*sizeof(double); // ode_args

    size += 1*nx*ns*sizeof(int); // ipiv

    size = (size + 63) / 64 * 64;
    size += 1 * 64;

    return size;
}

char *assign_irk_memory(sim_RK_opts *opts, sim_in *in, sim_irk_memory **memory, void *raw_memory)
{

    int nx = in->nx;
    int nu = in->nu; 
    // int NF = in->NF;
        
    int ns = opts->num_stages; // number of stages

    char *c_ptr = (char *)raw_memory;

    *memory = (sim_irk_memory *) c_ptr;
    c_ptr += sizeof(sim_irk_memory);

    (*memory)->JG = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->rG = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->K = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    struct d_strmat *JG = (*memory)->JG;
    struct d_strmat *rG = (*memory)->rG;
    struct d_strmat *K = (*memory)->K;

    // align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    d_create_strmat(nx*ns, nx*ns, JG, c_ptr);
    c_ptr += d_size_strmat(nx*ns, nx*ns);

    d_create_strmat(nx*ns, 1, rG, c_ptr);
    c_ptr += d_size_strmat(nx*ns, 1);

    d_create_strmat(nx*ns, 1, K, c_ptr);
    c_ptr += d_size_strmat(nx*ns, 1);

    (*memory)->xt0 = (double *)c_ptr;
    c_ptr += nx * sizeof(double);

    (*memory)->xt1 = (double *)c_ptr;
    c_ptr += nx * sizeof(double);

    (*memory)->kt = (double *)c_ptr;
    c_ptr += nx * sizeof(double);

    (*memory)->rGt = (double *)c_ptr;
    c_ptr += nx * sizeof(double);

    (*memory)->jac_out = (double *)c_ptr;
    c_ptr += 2 * nx * nx * sizeof(double);

    (*memory)->Jt = (double *)c_ptr;
    c_ptr += nx * nx * sizeof(double);

    (*memory)->ode_args = (double *)c_ptr;
    c_ptr += (2*nx+nu) * sizeof(double);

    (*memory)->ipiv = (int *)c_ptr;
    c_ptr += nx * ns * sizeof(int);

    return c_ptr;
}

sim_irk_memory *sim_irk_create_memory(sim_RK_opts *opts, sim_in *in)
{
    sim_irk_memory *memory;

    int_t bytes = irk_calculate_memory_size(opts, in);
    void *ptr = malloc(bytes);
    char *ptr_end = assign_irk_memory(opts, in, &memory, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    return memory;
}


int sim_irk_yt(const sim_in *in, sim_out *out, void *opts_, void *mem_){
    sim_RK_opts *opts = (sim_RK_opts *) opts_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    int ii, jj, iter, kk, ss;
    double a,b;

    int nx = in->nx;
	int nu = in->nu;
	// int NF = in->NF;
    int num_steps = in->num_steps;
    double *x = in->x;
	double *u = in->u;
    double step = in->step;
    // double *S_forw_in = in->S_forw;
    
    int num_stages = opts->num_stages;
    int newton_iter = opts->newton_iter;
    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
 
	double *xt0 = mem->xt0;
	double *xt1 = mem->xt1;
    double *kt = mem->kt;
    
	double *rGt = mem->rGt;
	double *jac_out = mem->jac_out;
	double *Jt = mem->Jt;
	int *ipiv = mem->ipiv;
	struct d_strmat *JG = mem->JG;
	struct d_strmat *rG = mem->rG;
    struct d_strmat *K = mem->K;

    double *xn = out->xn;

    printf("\nstep=%8.5f, stage=%d\n",step,num_stages);

    acados_timer timer;

    // initializae
    dgese_libstr(nx*num_stages, 1, 0.0, K, 0, 0);
    double *ode_args = mem->ode_args;
    for (ii=0;ii<nu;ii++)
        ode_args[2*nx+ii] = u[ii];

    struct d_strvec sxt0;
    d_create_strvec(nx, &sxt0, xt0);
    struct d_strvec sxt1;
    d_create_strvec(nx, &sxt1, xt1);
    struct d_strvec skt;
    d_create_strvec(nx, &skt, kt);

    // struct d_strvec sxt0, sxt1, skt;
    // void *mskt, *msxt0, *msxt1;
    // v_zeros_align(&msxt0, d_size_strvec(nx));
    // v_zeros_align(&msxt1, d_size_strvec(nx));
    // v_zeros_align(&mskt, d_size_strvec(nx));
    // d_create_strvec(nx, &sxt0, msxt0);
    // d_create_strvec(nx, &sxt1, msxt1);
    // d_create_strvec(nx, &skt, mskt);
    
    acados_tic(&timer);
    for(ss=0; ss<num_steps; ss++){

        for(iter=0; iter<newton_iter; iter++){
           
            for(ii=0; ii<num_stages; ii++){ // ii-th row of tableau   
                
                // d_cvt_vec2strvec(nx, x, &sxt0, 0);
                // d_cvt_vec2strvec(nx, x, &sxt1, 0);
                for(kk=0; kk<nx; kk++)
                    xt0[kk] = x[kk];

                // d_print_tran_strvec(nx, &sxt0, 0);

                for(jj=0; jj<num_stages; jj++){ // jj-th col of tableau
                    a = A_mat[ii+num_stages*jj];
                    // printf("%4.2f ",a);
					if(a!=0){
                        a *= step;
                        dcolex_libstr(nx, K, jj*nx, 0, &skt, 0);
                        // d_print_tran_strvec(nx, &skt, 0);
                        // daxpy_libstr(nx, a, K, jj*nx, &sxt0, 0, &sxt0, 0); // sxt0=sxt0+a*kj
                        daxpy_libstr(nx, a, &skt, 0, &sxt0, 0, &sxt0, 0); // sxt0=sxt0+a*kj
					}
                }
                // put xn+sum kj into first nx elements of ode_arg               
                d_cvt_strvec2vec(nx, &sxt0, 0, ode_args); 
                // put ki into the next nx elements of ode_arg

                // d_print_strmat(nx*num_stages, 1, K, 0, 0);                
                
                dcolex_libstr(nx, K, ii*nx, 0, &skt, 0);

                // d_print_tran_strvec(nx, &skt, 0);

                d_cvt_strvec2vec(nx, &skt, 0, ode_args+nx);      
                
                printf("\n");
                for (kk=0;kk<2*nx+nu;kk++)
                    printf("%5.3f ", ode_args[kk]);
                printf("\n");

                // compute the residual of implicit ode
                in->eval_impl_res(nx, nu, ode_args, rGt, in->impl_ode);   
                
                // printf("\n");
                // for (kk=0;kk<nx;kk++)
                //     printf("%5.3f ", rGt[kk]);
                // printf("\n");
                
                // compute the jacobian of implicit ode
                in->eval_impl_jac(nx, nu, ode_args, jac_out, in->impl_jac);

                // printf("\n");
                // // print_matrix("stdout",jac_out,nx,2*nx);
                // for(kk=0;kk<2*nx*nx;kk++)
                //     printf("%5.3f ", jac_out[kk]);
                // printf("\n");

                // compute the blocks of JG
                for (jj=0;jj<num_stages;jj++){
                    a = A_mat[ii+num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] = a* jac_out[kk];
                        if(jj==ii)
                            Jt[kk] += jac_out[nx*nx+kk];                          
                    }
                    // printf("\n");
                    // for (kk=0;kk< nx*nx;kk++)
                    //     printf("%5.3f ", Jt[kk]);
                    // printf("\n");

                    // fill in the block of JG
                    d_cvt_mat2strmat(nx, nx, Jt, nx, JG, ii*nx, jj*nx);
                    // fill in elements of rG
                    d_cvt_mat2strmat(nx, 1, rGt, nx, rG, ii*nx, 0);

                    // d_print_strmat(nx,nx,JG,ii*nx,jj*nx);
                } // end jj

            } // end ii
            
            // d_print_strmat(nx*num_stages, nx*num_stages, JG, 0, 0);
            // d_print_strmat(nx*num_stages, 1, rG, 0, 0);

            dgetrf_libstr(nx*num_stages, nx*num_stages, JG, 0, 0, JG, 0, 0, ipiv); // LU factorization with pivoting
            
            // printf("\n");
            // for (kk=0;kk<nx*num_stages;kk++)
            //     printf("%d ",ipiv[kk]);
            // printf("\n");

            // d_print_strmat(nx*num_stages,nx*num_stages,JG,0,0);

            drowpe_libstr(nx*num_stages, ipiv, rG);  // row permutations
            dtrsm_llnu_libstr(nx*num_stages, 1, 1.0, JG, 0, 0, rG, 0, 0, rG, 0, 0);  // L backsolve
            dtrsm_lunn_libstr(nx*num_stages, 1, 1.0, JG, 0, 0, rG, 0, 0, rG, 0, 0);  // U backsolve

            // d_print_strmat(nx*num_stages,1,rG,0,0);

            dgead_libstr(nx*num_stages, 1, -1.0, rG, 0, 0, K, 0, 0);

            // d_print_strmat(nx*num_stages,1, K, 0, 0);
        }// end iter

        for(ii=0;ii<num_stages;ii++){
            b = step * b_vec[ii];
            dcolex_libstr(nx, K, ii*nx, 0, &skt, 0);
            // daxpy_libstr(nx, b, &skt, 0, &sxt1, 0, &sxt1, 0);
            for(jj=0; jj<nx; jj++)
                x[jj] += b*kt[jj];
        }
        // d_cvt_strvec2vec(nx, &sxt1, 0, x);  

    }// end int step ss

    for (ii=0;ii<nx;ii++)
        xn[ii] = x[ii];

    out->info->CPUtime = acados_toc(&timer);

    // compute forward sens

    // v_free_align(msxt0);
    // v_free_align(msxt1);
    // v_free_align(mskt);
    // d_free_strvec(&sxt0);
    // d_free_strvec(&sxt1);
    // d_free_strvec(&skt);

    return 0;
}