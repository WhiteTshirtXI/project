/*  LBE: Lattice-Boltzmann Solver; Contains
		LBE_BCONDS (boundary conditions)
		LBE_MOVEXY (move velocity pointers for propagation in xy plane)
		LBE_ZCOL   (propagate in z direction and do collisions)  */
/***********************************************************************
Susp3D: Lattice-Boltzmann simulation code for particle-fluid suspensions
Copyright (C) 2003 Tony Ladd

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 ***********************************************************************/

#include "header.h"
extern int    max_x, max_y, max_z;
extern int    num_x, x_min, x_max;
extern int    num_proc, n_proc;
extern double tau[2], tau_v[2];

void lbe_bconds (Float ***velcs_df)
{
	extern int wall_flag;
	double send_buf[MAX_B], recv_buf[MAX_B];
	int x, y, z, xy, n, n_dir;
	int send_id, recv_id;
	int num_msg, num_buf, n_buf, buf_size, ptr;

	num_msg = MAX_B/(Num_Dir_X*(max_z+2));           /* # msgs per buffer */
	num_buf = ((max_y+2)-1)/num_msg + 1;             /* # buffers */

	send_id = (n_proc+1)%num_proc;               /* +x-direction (+1) */
	recv_id = (n_proc-1+num_proc)%num_proc;
	for (n_buf = 0; n_buf < num_buf; n_buf++)
	{
		buf_size = min((max_y+2)-n_buf*num_msg, num_msg);
		ptr  = 0;
		for (y = 0; y < buf_size; y++)
		{
			xy = num_x*(max_y+2) + y + n_buf*num_msg;
			for (n = 0; n < Num_Dir_X; n++)
			{
				n_dir = x_p[n];
				for (z = 0; z < (max_z+2); z++)
				{
					send_buf[ptr] = velcs_df[xy][n_dir][z];
					ptr++;
				}
			}
		}

		vector_xchg (send_buf, recv_buf, ptr, send_id, recv_id);           
		ptr  = 0;
		for (y = 0; y < buf_size; y++)
		{
			xy = y + n_buf*num_msg;
			for (n = 0; n < Num_Dir_X; n++)
			{
				n_dir = x_p[n];
				for (z = 0; z < max_z+2; z++)
				{
					velcs_df[xy][n_dir][z] = recv_buf[ptr];
					ptr++;
				}
			}
		}
	}

	send_id = (n_proc-1+num_proc)%num_proc;      /* -x-direction (-1) */
	recv_id = (n_proc+1)%num_proc;
	for (n_buf = 0; n_buf < num_buf; n_buf++)
	{
		buf_size = min(max_y+2-n_buf*num_msg, num_msg);
		ptr  = 0;
		for (y = 0; y < buf_size; y++)
		{
			xy = (max_y+2) + y + n_buf*num_msg;
			for (n = 0; n < Num_Dir_X; n++)
			{
				n_dir = x_m[n];
				for (z = 0; z < max_z+2; z++)
				{
					send_buf[ptr] = velcs_df[xy][n_dir][z];
					ptr++;
				}
			}
		}

		vector_xchg (send_buf, recv_buf, ptr, send_id, recv_id);           
		ptr  = 0;
		for (y = 0; y < buf_size; y++)
		{
			xy = num_x*(max_y+2) + (max_y+2) + y + n_buf*num_msg;
			for (n = 0; n < Num_Dir_X; n++)
			{
				n_dir = x_m[n];
				for (z = 0; z < max_z+2; z++)
				{
					velcs_df[xy][n_dir][z] = recv_buf[ptr];
					ptr++;
				}
			}
		}
	}

	/* Periodic boundary conditions for y and z*/
	/* +/- y direction */
	if(wall_flag < 1) 
		for(x=0; x<=num_x+1; x++) {
			xy = x*(max_y+2) + max_y; 
			for(n=0; n<Num_Dir_X; n++) {
				n_dir = y_p[n];
				for(z=0; z<=max_z+1; z++)
					velcs_df[x*(max_y+2)][n_dir][z] = velcs_df[xy][n_dir][z];
			}

			xy = x*(max_y+2) + 1;
			for(n=0; n<Num_Dir_X; n++) {
				n_dir = y_m[n];
				for(z=0; z<=max_z+1; z++)
					velcs_df[x*(max_y+2)+max_y+1][n_dir][z] = velcs_df[xy][n_dir][z];
			}
		}


	/* +/- z direction */
	if(wall_flag < 2)
		for(xy=0; xy<=(num_x+1)*(max_y+2)+max_y+1; xy++) {
			for(n=0; n<Num_Dir_X; n++) {
				n_dir = z_m[n];
				velcs_df[xy][n_dir][max_z+1]=velcs_df[xy][n_dir][1];

				n_dir = z_p[n];
				velcs_df[xy][n_dir][0]=velcs_df[xy][n_dir][max_z];
			}
		}
}


/*_MOVEXY: Updates velocity pointers for propagation in xy-plane
	Uses cyclic permutation of pointers */

void lbe_move (Float ***velcs_df)
{
	Float *tmp_p;
	Float ztmp, vztmp;
	int   x, y, z, xo, yo, q; 
	int   xy, xy_new, xy_old, z_new, z_old;
	int   n, n_dir;

	for (n = 0; n < Num_Dir_X; n++)              /* +x-direction */
		for (y = 0; y <= max_y+1; y++)
		{
			n_dir  = x_p[n];
			tmp_p  = velcs_df[num_x*(max_y+2)+y][n_dir];
			for (x = num_x; x >= 1; x--)
			{
				xy_new = x*(max_y+2) + y;
				xy_old = xy_new  - (max_y+2);
				velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
			}
			velcs_df[y][n_dir] = tmp_p;
		}

	for (n = 0; n < Num_Dir_X; n++)              /* -x-direction */
		for (y = 0; y <= max_y+1; y++)
		{
			n_dir  = x_m[n];
			tmp_p  = velcs_df[(max_y+2)+y][n_dir];
			for (x = 1; x <= num_x; x++)
			{
				xy_new = x*(max_y+2) + y;
				xy_old = xy_new  + (max_y+2);
				velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
			}
			velcs_df[(num_x+1)*(max_y+2)+y][n_dir] = tmp_p;
		}

	for (n = 0; n < Num_Dir_X; n++)              /* +y-direction */
		for (x = 0; x <= num_x+1; x++)
		{
			n_dir  = y_p[n];
			tmp_p  = velcs_df[x*(max_y+2)+max_y][n_dir];
			for (y = max_y; y >= 1; y--)
			{
				xy_new = x*(max_y+2) + y;
				xy_old = xy_new - 1;
				velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
			}
			velcs_df[x*(max_y+2)][n_dir] = tmp_p;
		}

	for (n = 0; n < Num_Dir_X; n++)              /* -y-direction */
		for (x = 0; x <= num_x+1; x++)
		{
			n_dir  = y_m[n];
			tmp_p  = velcs_df[x*(max_y+2)+1][n_dir];
			for (y = 1; y <= max_y; y++)
			{
				xy_new = x*(max_y+2) + y;
				xy_old = xy_new + 1;
				velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
			}
			velcs_df[x*(max_y+2)+(max_y+1)][n_dir] = tmp_p;
		}
	for (n = 0; n < Num_Dir_X; n++)               /* +z-direction */
		for (x=1; x<=num_x; x++) 
			for(y=1; y<=max_y; y++) {
				xy = x*(max_y+2)+y;

				n_dir  = z_p[n];
				ztmp  = velcs_df[xy][n_dir][max_z];
				for (z = max_z; z >= 1; z--) 
					velcs_df[xy][n_dir][z] = velcs_df[xy][n_dir][z-1];
				velcs_df[xy][n_dir][0] = ztmp;
			}

	for (n = 0; n < Num_Dir_X; n++) {              /* -z-direction */
		for (x=1; x<=num_x; x++) 
			for(y=1; y<=max_y; y++) {
				xy = x*(max_y+2)+y;
				n_dir  = z_m[n];
				ztmp  = velcs_df[xy][n_dir][1];
				for (z = 1; z <= max_z; z++) {
					vztmp = velcs_df[xy][n_dir][z+1];
					velcs_df[xy][n_dir][z] = velcs_df[xy][n_dir][z+1];
				}
				velcs_df[xy][n_dir][max_z+1] = ztmp;
			}
	}
}

/*  LBE_ZCOL: Collision operator and propagation in z for 1 column  */

void lbe_zcol(Float **velcz_df, int *node_map, struct vector f_ext, double sum_zcol[16], 
		int y, double rannum[], double fluidStress_pre[6], double fluidStress_pos[6])
{
	/*extern mt_state twister; */
	extern Float std_xi[Num_Dir], std_xi_1[Num_Dir];
	extern int add_noise;
	double modes[MAX_Z*10];
	double nfluc[Num_Dir], xi[Num_Dir];
	double m, m0, m1, m2, m3, m4, m5, m6, m7, m8, m9;
	double m1p, m2p, m3p, m4p, m5p, m6p, m7p, m8p, m9p;
	double m1m, m2m, m3m, m4m, m5m, m6m, m7m, m8m, m9m;
	double m0e, m4e, m5e, m6e, m7e, m8e, m9e;
	double lambda, lambda_v, lam[2], lam_v[2];
	int    z, z_p, z_m, z10, q;
	int    nyz;
	FILE   *stream;

	lam[0]   = 1.0 - 1.0/tau[0];                      /* Eigenvalues */
	lam_v[0] = 1.0 - 1.0/tau_v[0];
	lam[1]   = 1.0 - 1.0/tau[1];                      
	lam_v[1] = 1.0 - 1.0/tau_v[1];

	for (z10 = 0; z10 < 16; z10++)                 /* Zero sum_zcol array  */
		sum_zcol[z10] = 0.0;

	/*  Compute normal modes: m0 = mass; m1-m3 = momenta; m4-m9 = stress  */
	for (z = 1, z10 = 0; z <= max_z; z++, z10 += 10)
	{
		if (node_map[z] != 1)  {          /* Skip solid nodes */
			if(node_map[z] == 0) {
				lambda   = lam[0];                      /* Eigenvalues */
				lambda_v = lam_v[0];
				for(q=0; q<Num_Dir; q++)
					xi[q] = std_xi[q];
			}
			else {
				lambda   = lam[node_map[z]-1];                      /* Eigenvalues */
				lambda_v = lam_v[node_map[z]-1];
				for(q=0; q<Num_Dir; q++)
					xi[q] = std_xi_1[q];
			}

			nyz=((y-1)*max_z+(z-1))*6;

			m0  = velcz_df[ 0][z];
			m1p = velcz_df[ 1][z]   + velcz_df[ 2][z];
			m1m = velcz_df[ 1][z]   - velcz_df[ 2][z];
			m2p = velcz_df[ 3][z]   + velcz_df[ 4][z];
			m2m = velcz_df[ 3][z]   - velcz_df[ 4][z];
			m3p = velcz_df[ 5][z]   + velcz_df[ 6][z];
			m3m = velcz_df[ 5][z]   - velcz_df[ 6][z];
			m4p = velcz_df[ 7][z]   + velcz_df[ 8][z];
			m4m = velcz_df[ 7][z]   - velcz_df[ 8][z];
			m5p = velcz_df[ 9][z]   + velcz_df[10][z];
			m5m = velcz_df[ 9][z]   - velcz_df[10][z];
			m6p = velcz_df[11][z]   + velcz_df[12][z];
			m6m = velcz_df[11][z]   - velcz_df[12][z];
			m7p = velcz_df[13][z]   + velcz_df[14][z];
			m7m = velcz_df[13][z]   - velcz_df[14][z];
			m8p = velcz_df[15][z]   + velcz_df[16][z];
			m8m = velcz_df[15][z]   - velcz_df[16][z];
			m9p = velcz_df[17][z]   + velcz_df[18][z];
			m9m = velcz_df[17][z]   - velcz_df[18][z];

			m0 += m1p + m2p + m3p + m4p + m5p + m6p + m7p + m8p + m9p;
			m1  = m1m + m8m + m9m + m6m - m7m;
			m2  = m2m + m4m + m5m + m8m - m9m;
			m3  = m3m + m6m + m7m + m4m - m5m;
			m4  = m1p + m8p + m9p + m6p + m7p;
			m5  = m2p + m4p + m5p + m8p + m9p;
			m6  = m3p + m6p + m7p + m4p + m5p;
			m7  = m4p - m5p;
			m8  = m6p - m7p;
			m9  = m8p - m9p;

			m0e = m0 + Rho_Fl;                     /* Equilibrium distribution */
			m4e = m1*m1/m0e;
			m5e = m2*m2/m0e;
			m6e = m3*m3/m0e;
			m7e = m2*m3/m0e;
			m8e = m3*m1/m0e;
			m9e = m1*m2/m0e;
			m4 -= m4e + m0*CS2;                    /* Non-equilibrium stress */
			m5 -= m5e + m0*CS2;
			m6 -= m6e + m0*CS2;
			m7 -= m7e;
			m8 -= m8e;
			m9 -= m9e;

			sum_zcol[10] += m4/2.0;                /* Accumulate viscous stress */
			sum_zcol[11] += m5/2.0;
			sum_zcol[12] += m6/2.0;
			sum_zcol[13] += m7/2.0;
			sum_zcol[14] += m8/2.0;
			sum_zcol[15] += m9/2.0;

			fluidStress_pre[0] += m4/2.0;
			fluidStress_pre[1] += m5/2.0;
			fluidStress_pre[2] += m6/2.0;
			fluidStress_pre[3] += m9/2.0;
			fluidStress_pre[4] += m7/2.0;
			fluidStress_pre[5] += m8/2.0;

			m   = (m4 + m5 + m6)/3.0;              /* Collisions */
			m4 -= m;
			m5 -= m;
			m6 -= m;
			m4 *= lambda;
			m5 *= lambda;
			m6 *= lambda;
			m7 *= lambda;
			m8 *= lambda;
			m9 *= lambda;
			m  *= lambda_v;

			m4 += m;
			m5 += m;
			m6 += m;

			sum_zcol[10] += m4/2.0;                /* Accumulate viscous stress */
			sum_zcol[11] += m5/2.0;
			sum_zcol[12] += m6/2.0;
			sum_zcol[13] += m7/2.0;
			sum_zcol[14] += m8/2.0;
			sum_zcol[15] += m9/2.0;

			fluidStress_pos[0] += m4/2.0;
			fluidStress_pos[1] += m5/2.0;
			fluidStress_pos[2] += m6/2.0;
			fluidStress_pos[3] += m9/2.0;
			fluidStress_pos[4] += m7/2.0;
			fluidStress_pos[5] += m8/2.0;

			m4 += m4e;                             /* Non-equilibrium stress */
			m5 += m5e;
			m6 += m6e;
			m7 += m7e;
			m8 += m8e;
			m9 += m9e;

			/* Add extern stress */
			m4 +=  m1*f_ext.x*2.0/m0e;
			m5 +=  m2*f_ext.y*2.0/m0e;
			m6 +=  m3*f_ext.z*2.0/m0e;
			m7 += (m2*f_ext.z + m3*f_ext.y)/m0e;
			m8 += (m3*f_ext.x + m1*f_ext.z)/m0e;
			m9 += (m1*f_ext.y + m2*f_ext.x)/m0e;

			/* Add fluid fluctuations*/

			if(add_noise > 1) {
				m4 += xi[4]*rannum[nyz];
				m5 += xi[5]*rannum[nyz+1];
				m6 += xi[6]*rannum[nyz+2];
				m7 += xi[7]*rannum[nyz+3];
				m8 += xi[8]*rannum[nyz+4];
				m9 += xi[9]*rannum[nyz+5];
			}

			m1 += 0.5*f_ext.x;                     /* Add 1/2 external force */
			m2 += 0.5*f_ext.y;
			m3 += 0.5*f_ext.z;

			sum_zcol[0] += m0e;                    /* Accumulate mass and momenta */
			sum_zcol[1] += m1;
			sum_zcol[2] += m2;
			sum_zcol[3] += m3;
			sum_zcol[4] += m1*m1/m0e;              /* Accumulate Reynolds stress */
			sum_zcol[5] += m2*m2/m0e;
			sum_zcol[6] += m3*m3/m0e;
			sum_zcol[7] += m2*m3/m0e;
			sum_zcol[8] += m3*m1/m0e;
			sum_zcol[9] += m1*m2/m0e;

			m1 += 0.5*f_ext.x;                     /* Add 1/2 external force */
			m2 += 0.5*f_ext.y;
			m3 += 0.5*f_ext.z;

			modes[z10  ] = m0;                     /* Save modes for column */
			modes[z10+1] = m1;
			modes[z10+2] = m2;
			modes[z10+3] = m3;
			modes[z10+4] = m4;
			modes[z10+5] = m5;
			modes[z10+6] = m6;
			modes[z10+7] = m7;
			modes[z10+8] = m8;
			modes[z10+9] = m9;
		}
	}

	/*  Compute new population densities from modes  */

	for (z = 1, z10 = 0; z <= max_z; z++, z10 += 10)
	{
		if (node_map[z] != 1)            /* Fluid nodes */
		{
			m0  = modes[z10  ]/36.0;
			m1  = modes[z10+1]/12.0;
			m2  = modes[z10+2]/12.0;
			m3  = modes[z10+3]/12.0;
			m4  = modes[z10+4]/8.0;
			m5  = modes[z10+5]/8.0;
			m6  = modes[z10+6]/8.0;
			m7  = modes[z10+7]/4.0;
			m8  = modes[z10+8]/4.0;
			m9  = modes[z10+9]/4.0;

			m0 -=(m4 + m5 + m6)/3.0;
			m1p =(m0 + m4)*2.0;
			m2p =(m0 + m5)*2.0;
			m3p =(m0 + m6)*2.0;
			m4p = m0 + m5 + m6 + m7;
			m5p = m0 + m5 + m6 - m7;
			m6p = m0 + m6 + m4 + m8;
			m7p = m0 + m6 + m4 - m8;
			m8p = m0 + m4 + m5 + m9;
			m9p = m0 + m4 + m5 - m9;
			m1m = m1*2.0;
			m2m = m2*2.0;
			m3m = m3*2.0;
			m4m = m2 + m3;
			m5m = m2 - m3;
			m6m = m3 + m1;
			m7m = m3 - m1;
			m8m = m1 + m2;
			m9m = m1 - m2;

			velcz_df[ 0][z] = m0*12.0;
			velcz_df[ 1][z] = m1p + m1m;
			velcz_df[ 2][z] = m1p - m1m;
			velcz_df[ 3][z] = m2p + m2m;
			velcz_df[ 4][z] = m2p - m2m;
			velcz_df[ 5][z] = m3p + m3m;
			velcz_df[ 6][z] = m3p - m3m;
			velcz_df[ 7][z] = m4p + m4m;
			velcz_df[ 8][z] = m4p - m4m;
			velcz_df[ 9][z] = m5p + m5m;
			velcz_df[10][z] = m5p - m5m;
			velcz_df[11][z] = m6p + m6m;
			velcz_df[12][z] = m6p - m6m;
			velcz_df[13][z] = m7p + m7m;
			velcz_df[14][z] = m7p - m7m;
			velcz_df[15][z] = m8p + m8m;
			velcz_df[16][z] = m8p - m8m;
			velcz_df[17][z] = m9p + m9m;
			velcz_df[18][z] = m9p - m9m;
		}
		else                                       /* Zero solid nodes */
		{
			velcz_df[ 1][z] = 0.0;
			velcz_df[ 2][z] = 0.0;
			velcz_df[ 3][z] = 0.0;
			velcz_df[ 4][z] = 0.0;
			velcz_df[ 5][z] = 0.0;
			velcz_df[ 6][z] = 0.0;
			velcz_df[ 7][z] = 0.0;
			velcz_df[ 8][z] = 0.0;
			velcz_df[ 9][z] = 0.0;
			velcz_df[10][z] = 0.0;
			velcz_df[11][z] = 0.0;
			velcz_df[12][z] = 0.0;
			velcz_df[13][z] = 0.0;
			velcz_df[14][z] = 0.0;
			velcz_df[15][z] = 0.0;
			velcz_df[16][z] = 0.0;
			velcz_df[17][z] = 0.0;
			velcz_df[18][z] = 0.0;
		}
	}
}
//double* distribution_eq(int xy, int z, double ***velcs_df)
//{
//  double density = 0.;
//  struct vector vel = {0.,0.,0.};
//  //double density=Rho_Fl;
//  for(int q=0; q < 19; q++)
//  {
//    // density
//    density += velcs_df[xy][q][z];
//    // velocity
//    vel.x += velcs_df[xy][q][z]*c_x[q];
//    vel.y += velcs_df[xy][q][z]*c_y[q];
//    vel.z += velcs_df[xy][q][z]*c_z[q];
//  }
//  vel.x /= density;
//  vel.y /= density;
//  vel.z /= density;
//  double term3 = (vel.x*vel.x + vel.y*vel.y + vel.z*vel.z)*1.5;
//  double f_eq[19];
//  // distribution
//  for(int q=0; q < 19; q++)
//  {
//    //double term1 = (c_x[q]*velocity.x + c_y[q]*vel.y + c_z[q]*vel.z) / CS2;
//    double term1 = (c_x[q]*vel.x + c_y[q]*vel.y + c_z[q]*vel.z)*3.;
//    double term2 = 0.5*term1*term1;
//    f_eq[q] = weight[q]*density*(1 + term1 + term2 - term3); 
//  }
//  return f_eq;
//}

// Modification 20170410
void equilibrium_distrib(int xy, int z, double ***velcs_df, double dt,
     struct vector forceDen,  struct vector *correctedVel, double *f_eq)
{
  double density = 0.;
  correctedVel->x = 0.;
  correctedVel->y = 0.;
  correctedVel->z = 0.;
  
  for(int q=0; q < 19; q++) {
    density += velcs_df[xy][q][z];
    correctedVel->x += velcs_df[xy][q][z]*c_x[q];
    correctedVel->y += velcs_df[xy][q][z]*c_y[q];
    correctedVel->z += velcs_df[xy][q][z]*c_z[q];
  }
  // density += 36.0;
  //printf("density=%f\n",density);
  correctedVel->x += 0.5*dt*forceDen.x;
  correctedVel->y += 0.5*dt*forceDen.y;
  correctedVel->z += 0.5*dt*forceDen.z;
  correctedVel->x /= density;
  correctedVel->y /= density;
  correctedVel->z /= density;
    
  double term3 = (correctedVel->x*correctedVel->x + correctedVel->y*correctedVel->y + 
                  correctedVel->z*correctedVel->z)*1.5;
  // distribution
  for(int q=0; q < 19; q++)
  {
    double term1 = (c_x[q] * correctedVel->x + c_y[q] * correctedVel->y + 
                    c_z[q] * correctedVel->z)*3.;
    double term2 = 0.5*term1*term1;
    f_eq[q] = weight[q]*density*(1 + term1 + term2 - term3); 
  }  
}

void external_force(double tau, struct vector forceDen, struct vector correctedVel, 
     double *f_ext) 
{
//  double f_ext;
  double prefactor = 1 - 1. / (2*tau);
  double CS4 = CS2*CS2;
  for(int q = 0; q < 19; q++) {
    double term1 = c_x[q] * forceDen.x + c_y[q] * forceDen.y + c_z[q] * forceDen.z;
    double term2 = correctedVel.x*forceDen.x + correctedVel.y*forceDen.y + 
                   correctedVel.z*forceDen.z;
    double term3 = c_x[q]*correctedVel.x + c_y[q]*correctedVel.y + 
                   c_z[q]*correctedVel.z;
    f_ext[q] = prefactor*weight[q]*((term1-term2)/CS2 + (term3*term1/CS4));
  }
//  return f_ext;
}

// Modification 20170410
void collision(double tau, double ***velcs_df, struct vector **forceDen, double dt)
{
  int max_y_plus_2 = max_y+2;
  double tau_inverse = 1. / tau;
  for(int x=1; x <= max_x; x++) {
    for(int y=1; y <= max_y; y++) {
      int xy = x*(max_y_plus_2)+y;
      for(int z=1; z <= max_z; z++) {
        struct vector correctedVel;
        double f_eq[19];
        double f_ext[19];
        // because the structure of velcs_df, I can't use velcs_df[xy][z] as a parameter.
        equilibrium_distrib(xy, z, velcs_df, dt, forceDen[xy][z], &correctedVel, f_eq);
        external_force(tau, forceDen[xy][z], correctedVel, f_ext); 
        for(int q=0; q < 19; q++) {
          velcs_df[xy][q][z] = velcs_df[xy][q][z] + tau_inverse * 
                               (f_eq[q] - velcs_df[xy][q][z]) + f_ext[q]*dt;
        }
        // if(TURE)
        // fneq(x,y,z) u(x,y,z) forceDen(x,y,z)
        // for(int q=0; q<19; q++)
        // {
        //   +=fneq[xy][q][z]*c_y[q]*c_x[q];
        //   +=fneq[xy][q][z]*c_x[q]*c_y[q];
        //   +=fneq[xy][q][z]*c_x[q]*c_x[q];
        // }
        // *= -(1-0.5*dt/tau);
        // -=0.5*dt*(1-0.5*dt/tau)*(forceDen[xy][z].y*u.x+u.y*forceDen[xy][z].x);
      }
    }
  }
  // if(TRUE) 
  // volume average
}

void propagation(struct object *objects, int ** node_map, Float ***velcs_df)
{
  extern int num_obj;
  lbe_bconds(velcs_df);
  if(num_obj > 0)
    bnodes_wall(objects,velcs_df,node_map);
  lbe_move(velcs_df);
}

