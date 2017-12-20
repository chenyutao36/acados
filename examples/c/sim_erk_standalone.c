#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/crane_model/crane_model.h"

#define M_PI 3.14159265358979323846

int main() {

    int NX = 4;
    int NU = 1;
    int d = 1;

    double Ts = 0.05;
    sim_in sim_in;
    sim_out sim_out;
    sim_info info;
    sim_RK_opts rk_opts;
        
    sim_in.num_steps = 4;
    sim_in.step = Ts / sim_in.num_steps;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.sens_forw = true; // this has to be true
    sim_in.sens_adj = false;
    sim_in.sens_hess = false;
    sim_in.num_forw_sens = NX+NU;

    sim_in.vde = &vdeFun;
    sim_in.forward_vde_wrapper = &vde_fun;
    // sim_in.vde_adj = &adjFun;
    // sim_in.adjoint_vde_wrapper = &vde_adj_fun;

    sim_in.x = (double *)calloc( NX , sizeof(double) );
    sim_in.u = (double *)calloc( NU, sizeof(double) );
    sim_in.S_forw = (double *)calloc( NX * (NX + NU), sizeof(double) );
    sim_in.S_adj =(double *)calloc( NX + NU, sizeof(double) );  
    sim_in.grad_K = (double *)calloc(d * NX, sizeof(double) );

    sim_out.xn = (double *)calloc( NX , sizeof(double) );
    sim_out.S_forw = (double *)calloc( NX * (NX + NU), sizeof(double) );
    sim_out.info = &info;

    // initial state and control
    sim_in.x[1] = M_PI;
    for (int i=0;i <NU; i++) sim_in.u[i]=1.0;

    // initial sensitivity seed
    for (int i = 0; i < NX; i++)
        sim_in.S_forw[i * (NX + 1)] = 1.0;

    int workspace_size;
    sim_erk_create_arguments(&rk_opts, 4);
    workspace_size = sim_erk_calculate_workspace_size(&sim_in, &rk_opts);

    void *sim_work = (void *)malloc(workspace_size);

    sim_erk(&sim_in, &sim_out, &rk_opts, 0, sim_work);

    double *xn = sim_out.xn;

    printf("\nxn: \n");
    for (int ii=0;ii<NX;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    // if you don't need sensitivities, just ignore 
    // double *S_forw_out = sim_out.S_forw;  
    // if (sim_in.sens_forw){
    //     printf("\nS_forw_out: \n");
    //     for (int ii=0;ii<NX;ii++){
    //         for (int jj=0;jj<NX+NU;jj++)
    //             printf("%8.5f ",S_forw_out[jj*NX+ii]);
    //         printf("\n");
    //     }
    // }

}