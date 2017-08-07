/*	CLUSTER_FORCE: Subtracts singular lubrication forces from clusters */
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

void cluster_force (struct object *objects, struct cluster *clusters, double dt, int num_cl, int sign)
{
	extern double tau;
	extern int 	max_x, max_y, max_z;
	extern int 	num_sph, num_obj;
	double      x12, y12, z12, r12;
    double      Fx, Fy, Fz;
    double      S1xx, S1yy, S1zz, S1yz, S1zx, S1xy;
    double      S2xx, S2yy, S2zz, S2yz, S2zx, S2xy;
    double      Xa, Xa_max, Xg1, Xg2;
    double      nu, rad, dr, cut;
    double      alpha, umx, umy, umz, um_r;
	double      beta1, beta2, beta3, beta4, beta5; 
    int         n1_sph, n2_obj;
	int         n1, n2, n_cl;

    nu = (tau - 0.5)/3.0;

	for (n_cl = 0; n_cl < num_cl; n_cl++)
	for (n1 = 0; n1 < clusters[n_cl].list[0]; n1++)
    {
		n1_sph = clusters[n_cl].list[n1+1];
		for (n2 = 0; n2 < objects[n1_sph].cl_list[0]; n2++)
       	{
           	n2_obj  = objects[n1_sph].cl_list[n2+1];
			if (n2_obj < n1_sph)  continue;                                                 /* Skip pairs with n2 < n1 */

			if (n2_obj < num_sph)
			{
				x12 = n_image(objects[n1_sph].r.x - objects[n2_obj].r.x, max_x);
				y12 = n_image(objects[n1_sph].r.y - objects[n2_obj].r.y, max_y);
				z12 = n_image(objects[n1_sph].r.z - objects[n2_obj].r.z, max_z);
				r12 = sqrt(x12*x12 + y12*y12 + z12*z12);
				rad = (objects[n1_sph].r_a + objects[n2_obj].r_a)/2.0;
				dr  = r12 - 2.0*rad;
				umx = objects[n1_sph].u.x;
				umy = objects[n1_sph].u.y;
				umz = objects[n1_sph].u.z;
				if (objects[n2_obj].mass_flag)
				{
					umx -= objects[n2_obj].u.x;                                             /* Include velocity of finite mass particles */
					umy -= objects[n2_obj].u.y;
					umz -= objects[n2_obj].u.z;
				}
			}
			else
			{
				x12 = (objects[n1_sph].r.x - objects[n2_obj].r.x)*objects[n2_obj].e.x*objects[n2_obj].e.x;
				y12 = (objects[n1_sph].r.y - objects[n2_obj].r.y)*objects[n2_obj].e.y*objects[n2_obj].e.y;
				z12 = (objects[n1_sph].r.z - objects[n2_obj].r.z)*objects[n2_obj].e.z*objects[n2_obj].e.z;
				r12 = sqrt(x12*x12 + y12*y12 + z12*z12);
				rad = objects[n1_sph].r_a;
				dr  = r12 - rad;
				umx = objects[n1_sph].u.x;
				umy = objects[n1_sph].u.y;
				umz = objects[n1_sph].u.z;
			} 

           	alpha = dt*Pi*Rho_Fl*nu*rad;
			cut = (objects[n1_sph].lub_cut + objects[n2_obj].lub_cut)/2.0;
			dr  = max(dr, cut);
			if (n2_obj < num_sph)
			{
				beta1 = objects[n2_obj].r_a/objects[n1_sph].r_a;
				beta2 = beta1*beta1;
				beta3 = beta1*beta2;
				beta4 = pow((1 + beta1), 4.0);
				beta5 = pow((1 + beta1), 5.0);
				Xa    =  4.0*beta2/beta4;
				Xg1   = 12.0*beta2/beta5;
				Xg2   =-12.0*beta3/beta5;
			}
			else
			{
				Xa    =  1.0;
				Xg1   =  3.0/2.0;
				Xg2   =  0.0;
			}
          	Xa  *= rad/dr;
			Xg1 *= rad/dr;
			Xg2 *= rad/dr;
          	Xa  *= 6.0*alpha;
			Xg1 *= 4.0*alpha*rad;
			Xg2 *= 4.0*alpha*rad;

			if (objects[n2_obj].mass_flag)  Xa_max = Kappa*min(objects[n1_sph].mass, objects[n2_obj].mass);
			else                            Xa_max = Kappa*objects[n1_sph].mass;
			if (n2_obj < num_sph)
			{
				Xa  -= Xa_max;
				Xg1 -= Xa_max*rad*2.0/(1.0+beta1);
				Xg2 += Xa_max*rad*2.0*beta1/(1.0+beta1);
			}
			else
			{
				Xa  -= Xa_max;
				Xg1 -= Xa_max*rad;
			}
			
			x12  /= r12;
			y12  /= r12;
			z12  /= r12;
			um_r  = umx*x12 + umy*y12 + umz*z12;

           	Fx    = -Xa*um_r*x12;
           	Fy    = -Xa*um_r*y12;
          	Fz    = -Xa*um_r*z12;
			S1xx  =  Xg1*um_r*x12*x12;
			S1yy  =  Xg1*um_r*y12*y12;
			S1zz  =  Xg1*um_r*z12*z12;
			S1yz  =  Xg1*um_r*y12*z12;
			S1zx  =  Xg1*um_r*z12*x12;
			S1xy  =  Xg1*um_r*x12*y12;
			S2xx  = -Xg2*um_r*x12*x12;
			S2yy  = -Xg2*um_r*y12*y12;
			S2zz  = -Xg2*um_r*z12*z12;
			S2yz  = -Xg2*um_r*y12*z12;
			S2zx  = -Xg2*um_r*z12*x12;
			S2xy  = -Xg2*um_r*x12*y12;

           	objects[n1_sph].f.x += sign*Fx;
           	objects[n1_sph].f.y += sign*Fy;
           	objects[n1_sph].f.z += sign*Fz;
           	objects[n2_obj].f.x -= sign*Fx;
           	objects[n2_obj].f.y -= sign*Fy;
           	objects[n2_obj].f.z -= sign*Fz;
			objects[n1_sph].pf_lub.xx += sign*S1xx;
			objects[n1_sph].pf_lub.yy += sign*S1yy;
			objects[n1_sph].pf_lub.zz += sign*S1zz;
			objects[n1_sph].pf_lub.yz += sign*S1yz;
			objects[n1_sph].pf_lub.zx += sign*S1zx;
			objects[n1_sph].pf_lub.xy += sign*S1xy;
			objects[n2_obj].pf_lub.xx += sign*S2xx;
			objects[n2_obj].pf_lub.yy += sign*S2yy;
			objects[n2_obj].pf_lub.zz += sign*S2zz;
			objects[n2_obj].pf_lub.yz += sign*S2yz;
			objects[n2_obj].pf_lub.zx += sign*S2zx;
			objects[n2_obj].pf_lub.xy += sign*S2xy;
		}
   	}
}
