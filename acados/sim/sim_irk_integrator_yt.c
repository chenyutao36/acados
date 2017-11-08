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
int NF = in->NF;
 
int ns = opts->num_stages; // number of stages

int size = sizeof(sim_irk_memory);

size += 1*sizeof(struct d_strmat); // JG
size += 2*sizeof(struct d_strvec); // rG, K

size += 1*d_size_strmat(nx*ns, nx*ns); // JG
size += 2*d_size_strvec(nx*ns); // rG K

size += 3*nx*sizeof(double); // x, xt0, xt1
size += 1*nu*sizeof(double); // u
size += 2*nx*ns*sizeof(double); // Kt, rGt
size += 2*nx*(2*nx+nu)*sizeof(double); // Jt0 Jt1

size += 1*nx*ns*sizeof(int); // ipiv

size = (size + 63) / 64 * 64;
size += 1 * 64;

return size;
}

int sim_irk_yt(const sim_in *in, sim_out *out, void *opts_, void *mem_){
    sim_RK_opts *opts = (sim_RK_opts *) opts_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    int ii, jj, step, iter, ss;
    double a,b;

    int nx = in->nx;
	int nu = in->nu;
	int NF = in->NF;
    int num_steps = in->num_steps;
    double step = in->step;
    double *S_forw_in = in->S_forw;
    
    int num_stages = opts->num_stages;
    int newton_iter = opts->newton_iter;
    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;

    double *x = mem->x;
	double *u = mem->u;
	double *xt0 = mem->xt0;
	double *xt1 = mem->xt1;
	double *Kt = mem->Kt;
	double *rGt = mem->rGt;
	double *Jt0 = mem->Jt0;
	double *Jt1 = mem->Jt1;
	int *ipiv = mem->ipiv;
	struct d_strmat *JG = mem->JG;
	struct d_strvec *rG = mem->rG;
    struct d_strvec *K = mem->K;

    // test
    dvecse_libstr(nx*num_stages, 0.0, K, 0);
    double *ode_args = mem->ode_args;
    for (ii=0;ii<nu;ii++)
        ode_args[2*nx+ii] = u[ii];

    struct d_strvec sxt0;
	d_create_strvec(nx, &sxt0, xt0);
	struct d_strvec sxt1;
    d_create_strvec(nx, &sxt1, xt1);
    

    for(step=0; step<num_steps; step++){

        for(iter=0; iter<newton_iter; iter++){

            for(ii=0; ii<num_stages; ss++){

                // extract initial states into xt0
                for(jj=0; jj<nx; jj++)
                    xt0[jj] = x[jj];

                for(jj=0; jj<num_stages; jj++){
					a = A_mat[ss+num_stages*jj];
					if(a!=0){
						a *= step;
						daxpy_libstr(nx, a, &K, jj*nx, &sxt0, 0, &sxt0, 0); // sxt0 = sxt0 + a * sxt1
					}
                }
                // put xn+sum kj into first nx elements of ode_arg               
                d_cvt_strvec2vec(nx, &sxt0, 0, ode_args); 
                // put ki into the next nx elements of ode_arg
                d_cvt_strvec2vec(nx, &K, ss*nx, ode_args+nx);            
                // compute the residual of implicit ode
                in->impl_ode_fun(nx,nu,ode_args,rGt+ii*nx,in->impl_ode);
                // compute the jacobian of implicit ode
                in->impl_ode_jac_fun(nx,nu,ode_args,Jt0,in->impl_ode_jac);

                for (jj=0;jj<num_stages;jj++){
                    a = A_mat[ss+num_stages*jj];
                    if (a!=0){
                        a *= step;
                    }
                    
                }
            }


        }

    }
}