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

#ifndef ACADOS_SIM_SIM_ERK_INTEGRATOR_YT_H_
#define ACADOS_SIM_SIM_ERK_INTEGRATOR_YT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_rk_common_new.h"
#include "acados/sim/sim_common_new.h"
#include "acados/utils/types.h"

typedef struct{  

    real_t *rhs_forw_in;  // x + S + p

    real_t *K_traj; // (stages *nX) or (steps*stages*nX) for adj
    real_t *out_forw_traj; // S or (steps+1)*nX for adj

    real_t *rhs_adj_in;
    real_t *out_adj_tmp;
    real_t *adj_traj;    
    
} sim_erk_memory;

int_t erk_calculate_memory_size(sim_RK_opts *opts, sim_in *in);

char *assign_erk_memory(sim_RK_opts *opts, sim_in *in, sim_erk_memory **memory, void *raw_memory);

sim_erk_memory *sim_erk_create_memory(sim_RK_opts *opts, sim_in *in);

int_t sim_erk_new(const sim_in *in, sim_out *out, void *opts_, void *mem_);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_