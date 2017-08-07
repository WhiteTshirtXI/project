/*  LUB: Calculates lubrication forces and stresses for polydisperse spheres */
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

void lub  (struct object *objects, double dt)
{
	extern double tau;
	extern double lub_N, lub_T, lub_R;
	extern int    max_x, max_y, max_z;
	extern int    num_sph, num_obj;
	double        x12, y12, z12, r12;
	double        umx, umy, umz, w1x, w1y, w1z, w2x, w2y, w2z;
	double        umxrx, umxry, umxrz, w1xrx, w1xry, w1xrz, w2xrx, w2xry, w2xrz;
	double        um_r, w1_r, w2_r;
	double        Xa, Ya, Yb1, Yb2, Xg1, Xg2, Yg1, Yg2;
	double        Yc11, Yc12, Yc21, Yc22, Yh11, Yh12, Yh21, Yh22;
	double        Fx, Fy, Fz, T1x, T1y, T1z, T2x, T2y, T2z;
	double        S1xx, S1yy, S1zz, S1yz, S1zx, S1xy;
	double        S2xx, S2yy, S2zz, S2yz, S2zx, S2xy;
	double        nu, rad, cut, dr, h, h_c;
	double        alpha, beta, psi;
	double        beta1, beta2, beta3, beta4, beta5;
	int           n1_sph, n2_sph, n2_obj, n;


	for (n1_sph = 0; n1_sph < num_sph; n1_sph++)		               /* Initialize cluster lists */
		 objects[n1_sph].cl_list[0] = 0;

	for (n1_sph = 1; n1_sph < num_sph; n1_sph++)                       /* Spheroid-Spheroid lubrication */
	{	
		for (n = 1; n <= objects[n1_sph].list[0]; n++)
		{
			n2_sph  = objects[n1_sph].list[n];			
			if (n2_sph < n1_sph)
			{
				cut = (objects[n1_sph].lub_cut + objects[n2_sph].lub_cut)/2.0;
				rad = (objects[n1_sph].r_a + objects[n2_sph].r_a)/2.0;
				x12 = n_image(objects[n1_sph].r.x - objects[n2_sph].r.x, max_x);
				y12 = n_image(objects[n1_sph].r.y - objects[n2_sph].r.y, max_y);
				z12 = n_image(objects[n1_sph].r.z - objects[n2_sph].r.z, max_z);
				r12 = sqrt(x12*x12 + y12*y12 + z12*z12);
				dr  = r12-2.0*rad;
				h   = max(dr, cut);

				if (dr < max(cut, max(lub_N, max(lub_T, lub_R))))
				{
					nu    = (tau - 0.5)/3.0;
					alpha = Pi*Rho_Fl*nu*rad*dt;
					beta1 = objects[n2_sph].r_a/objects[n1_sph].r_a;
					beta2 = beta1*beta1;
					beta3 = beta1*beta2;
					beta4 = pow((1 + beta1), 4.0);
					beta5 = pow((1 + beta1), 5.0);
					if (dr < max(cut, lub_N))
					{
						Xa   =  4.0*beta2/beta4;                         /* Set lubrication coefficients */
						Xg1  = 12.0*beta2/beta5;
						Xg2  =-12.0*beta3/beta5;
						h_c   = max(dr, lub_N);
						Xa   *= rad/h - rad/h_c;
						Xg1  *= rad/h - rad/h_c;
						Xg2  *= rad/h - rad/h_c;
						Xa   *= 6.0*alpha;
						Xg1  *= 4.0*alpha*rad;
						Xg2  *= 4.0*alpha*rad;
					}
					else
						Xa = Xg1 = Xg2 = 0.0;
					if (dr < max(cut, lub_T))
					{
						Ya    = 8.0*beta1*(2.0+beta1+2.0*beta2)/(15.0*beta4);
						Yb1   =-4.0*beta1*(4.0+beta1)/(5.0*beta4);
						Yb2   = 4.0*beta1*(4.0*beta2+beta1)/(5.0*beta4);
						Yg1   = 2.0*beta1*(4.0-beta1+7.0*beta2)/(5.0*beta5);
						Yg2   =-2.0*beta2*(4.0*beta2-beta1+7.0)/(5.0*beta5);
						h_c   = max(dr, lub_T);
						beta  = alpha*log(h_c/h);
						Ya   *= 6.0*beta;
						Yb1  *= 4.0*beta*rad;
						Yb2  *= 4.0*beta*rad;
						Yg1  *= 4.0*beta*rad;
						Yg2  *= 4.0*beta*rad;
					}
					else
						Ya = Yb1 = Yb2 = Yg1 = Yg2 = 0.0;
					if (dr < max(cut, lub_R))
					{
						Yc11  = 16.0*beta1/(5.0*beta4);
						Yc12  =  4.0*beta2/(5.0*beta4);
						Yc21  =  4.0*beta2/(5.0*beta4);
						Yc22  = 16.0*beta3/(5.0*beta4);
						Yh11  =  4.0*beta1*(2.0-beta1)/(5.0*beta5);
						Yh12  =  2.0*beta2*(1.0+7.0*beta1)/(5.0*beta5);
						Yh21  =  2.0*beta2*(beta1+7.0)/(5.0*beta5);
						Yh22  =  4.0*beta2*(2.0*beta2-beta1)/(5.0*beta5);
						h_c   = max(dr, lub_R);
						beta  = alpha*log(h_c/h);
						Yc11 *= 8.0*beta*rad*rad;
						Yc12 *= 8.0*beta*rad*rad;
						Yc21 *= 8.0*beta*rad*rad;
						Yc22 *= 8.0*beta*rad*rad;
						Yh11 *= 8.0*beta*rad*rad;
						Yh21 *= 8.0*beta*rad*rad;
						Yh12 *= 8.0*beta*rad*rad;
						Yh22 *= 8.0*beta*rad*rad;
					}
					else
						Yc11 = Yc12 = Yc21 = Yc22 = Yh11 = Yh12 = Yh21 = Yh22 = 0.0;

					x12  /= r12;
					y12  /= r12;
					z12  /= r12;
					umx   = objects[n1_sph].u.x - objects[n2_sph].u.x;
					umy   = objects[n1_sph].u.y - objects[n2_sph].u.y;
					umz   = objects[n1_sph].u.z - objects[n2_sph].u.z;
					w1x   = objects[n1_sph].w.x;
					w1y   = objects[n1_sph].w.y;
					w1z   = objects[n1_sph].w.z;
					w2x   = objects[n2_sph].w.x;
					w2y   = objects[n2_sph].w.y;
					w2z   = objects[n2_sph].w.z;
					um_r  = umx*x12 + umy*y12 + umz*z12;
					w1_r  = w1x*x12 + w1y*y12 + w1z*z12;
					w2_r  = w2x*x12 + w2y*y12 + w2z*z12;
					umxrx = umy*z12 - umz*y12;
					umxry = umz*x12 - umx*z12;
					umxrz = umx*y12 - umy*x12;
					w1xrx = w1y*z12 - w1z*y12;
					w1xry = w1z*x12 - w1x*z12;
					w1xrz = w1x*y12 - w1y*x12;
					w2xrx = w2y*z12 - w2z*y12;
					w2xry = w2z*x12 - w2x*z12;
					w2xrz = w2x*y12 - w2y*x12;
					
					Fx    = -Xa*um_r*x12 - Ya*(umx-um_r*x12);
					Fy    = -Xa*um_r*y12 - Ya*(umy-um_r*y12);
					Fz    = -Xa*um_r*z12 - Ya*(umz-um_r*z12);
					Fx   += -Yb1*w1xrx + Yb2*w2xrx;
					Fy   += -Yb1*w1xry + Yb2*w2xry;
					Fz   += -Yb1*w1xrz + Yb2*w2xrz;
					T1x   =  Yb1*umxrx;
					T1y   =  Yb1*umxry;
					T1z   =  Yb1*umxrz;
					T2x   = -Yb2*umxrx;
					T2y   = -Yb2*umxry;
					T2z   = -Yb2*umxrz;
					T1x  += -Yc11*(w1x-w1_r*x12) - Yc12*(w2x-w2_r*x12);
					T1y  += -Yc11*(w1y-w1_r*y12) - Yc12*(w2y-w2_r*y12);
					T1z  += -Yc11*(w1z-w1_r*z12) - Yc12*(w2z-w2_r*z12);
					T2x  += -Yc21*(w1x-w1_r*x12) - Yc22*(w2x-w2_r*x12);
					T2y  += -Yc21*(w1y-w1_r*y12) - Yc22*(w2y-w2_r*y12);
					T2z  += -Yc21*(w1z-w1_r*z12) - Yc22*(w2z-w2_r*z12);
					S1xx  =  Xg1*um_r*x12*x12 + Yg1*(x12*umx+umx*x12);
					S1yy  =  Xg1*um_r*y12*y12 + Yg1*(y12*umy+umy*y12);
					S1zz  =  Xg1*um_r*z12*z12 + Yg1*(z12*umz+umz*z12);
					S1yz  =  Xg1*um_r*y12*z12 + Yg1*(y12*umz+umy*z12);
					S1zx  =  Xg1*um_r*z12*x12 + Yg1*(z12*umx+umz*x12);
					S1xy  =  Xg1*um_r*x12*y12 + Yg1*(x12*umy+umx*y12);
					S2xx  = -Xg2*um_r*x12*x12 - Yg2*(x12*umx+umx*x12);
					S2yy  = -Xg2*um_r*y12*y12 - Yg2*(y12*umy+umy*y12);
					S2zz  = -Xg2*um_r*z12*z12 - Yg2*(z12*umz+umz*z12);
					S2yz  = -Xg2*um_r*y12*z12 - Yg2*(y12*umz+umy*z12);
					S2zx  = -Xg2*um_r*z12*x12 - Yg2*(z12*umx+umz*x12);
					S2xy  = -Xg2*um_r*x12*y12 - Yg2*(x12*umy+umx*y12);
					S1xx += -Yh11*(x12*w1xrx+w1xrx*x12) - Yh12*(x12*w2xrx+w2xrx*x12);
					S1yy += -Yh11*(y12*w1xry+w1xry*y12) - Yh12*(y12*w2xry+w2xry*y12);
					S1zz += -Yh11*(z12*w1xrz+w1xrz*z12) - Yh12*(z12*w2xrz+w2xrz*z12);
					S1yz += -Yh11*(y12*w1xrz+w1xry*z12) - Yh12*(y12*w2xrz+w2xry*z12);
					S1zx += -Yh11*(z12*w1xrx+w1xrz*x12) - Yh12*(z12*w2xrx+w2xrz*x12);
					S1xy += -Yh11*(x12*w1xry+w1xrx*y12) - Yh12*(x12*w2xry+w2xrx*y12);
					S2xx += -Yh21*(x12*w1xrx+w1xrx*x12) - Yh22*(x12*w2xrx+w2xrx*x12);
					S2yy += -Yh21*(y12*w1xry+w1xry*y12) - Yh22*(y12*w2xry+w2xry*y12);
					S2zz += -Yh21*(z12*w1xrz+w1xrz*z12) - Yh22*(z12*w2xrz+w2xrz*z12);
					S2yz += -Yh21*(y12*w1xrz+w1xry*z12) - Yh22*(y12*w2xrz+w2xry*z12);
					S2zx += -Yh21*(z12*w1xrx+w1xrz*x12) - Yh22*(z12*w2xrx+w2xrz*x12);
					S2xy += -Yh21*(x12*w1xry+w1xrx*y12) - Yh22*(x12*w2xry+w2xrx*y12);

					objects[n1_sph].f.x += Fx;			
					objects[n1_sph].f.y += Fy;
					objects[n1_sph].f.z += Fz;
					objects[n2_sph].f.x -= Fx;			
					objects[n2_sph].f.y -= Fy;
					objects[n2_sph].f.z -= Fz;
					objects[n1_sph].t.x += T1x;
					objects[n1_sph].t.y += T1y;
					objects[n1_sph].t.z += T1z;
					objects[n2_sph].t.x += T2x;
					objects[n2_sph].t.y += T2y;
					objects[n2_sph].t.z += T2z;
					objects[n1_sph].pf_lub.xx += S1xx;
					objects[n1_sph].pf_lub.yy += S1yy;
					objects[n1_sph].pf_lub.zz += S1zz;
					objects[n1_sph].pf_lub.yz += S1yz;
					objects[n1_sph].pf_lub.zx += S1zx;
					objects[n1_sph].pf_lub.xy += S1xy;
					objects[n2_sph].pf_lub.xx += S2xx;
					objects[n2_sph].pf_lub.yy += S2yy;
					objects[n2_sph].pf_lub.zz += S2zz;
					objects[n2_sph].pf_lub.yz += S2yz;
					objects[n2_sph].pf_lub.zx += S2zx;
					objects[n2_sph].pf_lub.xy += S2xy;


					switch (objects[n1_sph].mass_flag + objects[n2_sph].mass_flag*2)
					{
						case 0:
							psi = 0.0;
						break;
						case 1:
							psi = Xa/objects[n1_sph].mass;
						break;
						case 2:
							psi = Xa/objects[n2_sph].mass;
						break;
						case 3:
							psi = Xa/min(objects[n1_sph].mass, objects[n2_sph].mass);
						break;
					}
					if (psi > Kappa)                                   /* List pairs for implicit update */
					{
						objects[n1_sph].cl_list[0]++;
						objects[n2_sph].cl_list[0]++;
						objects[n1_sph].cl_list[objects[n1_sph].cl_list[0]] = n2_sph;
						objects[n2_sph].cl_list[objects[n2_sph].cl_list[0]] = n1_sph;
						if (objects[n1_sph].cl_list[0] > MAX_N ) fatal_err("too many neighbors in cluster", MAX_N);
						if (objects[n2_sph].cl_list[0] > MAX_N ) fatal_err("too many neighbors in cluster", MAX_N);
					}
				}
			}
		}
	}

	for (n1_sph = 0; n1_sph <  num_sph; n1_sph++)                          /* Spheroid-Wall lubrication */
	{
		for (n2_obj = num_sph; n2_obj < num_obj; n2_obj++)
		{
			rad =  objects[n1_sph].r_a;
			cut = (objects[n1_sph].lub_cut + objects[n2_obj].lub_cut)/2.0;
			x12 = (objects[n1_sph].r.x - objects[n2_obj].r.x)*objects[n2_obj].e.x*objects[n2_obj].e.x;
			y12 = (objects[n1_sph].r.y - objects[n2_obj].r.y)*objects[n2_obj].e.y*objects[n2_obj].e.y;
			z12 = (objects[n1_sph].r.z - objects[n2_obj].r.z)*objects[n2_obj].e.z*objects[n2_obj].e.z;
			r12 = sqrt(x12*x12 + y12*y12 + z12*z12);
			dr  = r12-rad;
			h   = max(dr, cut);

			if (dr < max(cut, max(lub_N, max(lub_T, lub_R))))
			{
				nu    = (tau - 0.5)/3.0;
				alpha = Pi*Rho_Fl*nu*rad*dt;
				if (dr < max(cut, lub_N))
				{
					Xa    = 1.0;                                             /* Set lubrication coefficients */
					Xg1   = 3.0/2.0;
					h_c   = max(dr, lub_N);
					Xa   *= rad/h - rad/h_c;
					Xg1  *= rad/h - rad/h_c;
					Xa   *= 6.0*alpha;
					Xg1  *= 4.0*alpha*rad;
				}
				else
					Xa = Xg1 = 0.0;
				if (dr < max(cut, lub_T))
				{
					Ya    = 8.0/15.0;
					Yb1   =-1.0/5.0;
					Yg1   = 7.0/10.0;
					h_c   = max(dr, lub_T);
					beta  = alpha*log(h_c/h);
					Ya   *= 6.0*beta;
					Yb1  *= 4.0*beta*rad;
					Yg1  *= 4.0*beta*rad;
				}
				else
					Ya = Yb1 = Yg1 = 0.0;
				if (dr < max(cut, lub_R))
				{
					Yc11  = 2.0/5.0;
					Yh11  =-1.0/10.0;
					h_c   = max(dr, lub_R);
					beta  = alpha*log(h_c/h);
					Yc11 *= 8.0*beta*rad*rad;
					Yh11 *= 8.0*beta*rad*rad;
				}
				else 
					Yc11 = Yh11 = 0.0;

				x12 /= r12;
				y12 /= r12;
				z12 /= r12;
				umx = objects[n1_sph].u.x - objects[n2_obj].u.x;
				umy = objects[n1_sph].u.y - objects[n2_obj].u.y;
				umz = objects[n1_sph].u.z - objects[n2_obj].u.z;
				w1x = objects[n1_sph].w.x;
				w1y = objects[n1_sph].w.y;
				w1z = objects[n1_sph].w.z;
				um_r = umx*x12 + umy*y12 + umz*z12;
				w1_r = w1x*x12 + w1y*y12 + w1z*z12;
				umxrx = umy*z12 - umz*y12;
				umxry = umz*x12 - umx*z12;
				umxrz = umx*y12 - umy*x12;
				w1xrx = w1y*z12 - w1z*y12;
				w1xry = w1z*x12 - w1x*z12;
				w1xrz = w1x*y12 - w1y*x12;
					
				Fx    = -Xa*x12*um_r - Ya*(umx-um_r*x12);
				Fy    = -Xa*y12*um_r - Ya*(umy-um_r*y12);
				Fz    = -Xa*z12*um_r - Ya*(umz-um_r*z12);
				Fx   += -Yb1*w1xrx;
				Fy   += -Yb1*w1xry;
				Fz   += -Yb1*w1xrz;
				T1x   =  Yb1*umxrx;
				T1y   =  Yb1*umxry;
				T1z   =  Yb1*umxrz;
				T1x  += -Yc11*(w1x-w1_r*x12);
				T1y  += -Yc11*(w1y-w1_r*y12);
				T1z  += -Yc11*(w1z-w1_r*z12);
				S1xx  =  Xg1*(um_r*x12*x12) + Yg1*(x12*umx+umx*x12);
				S1yy  =  Xg1*(um_r*y12*y12) + Yg1*(y12*umy+umy*y12);
				S1zz  =  Xg1*(um_r*z12*z12) + Yg1*(z12*umz+umz*z12);
				S1yz  =  Xg1*(um_r*y12*z12) + Yg1*(y12*umz+umy*z12);
				S1zx  =  Xg1*(um_r*z12*x12) + Yg1*(z12*umx+umz*x12);
				S1xy  =  Xg1*(um_r*x12*y12) + Yg1*(x12*umy+umx*y12);
				S1xx += -Yh11*(x12*w1xrx+w1xrx*x12);
				S1yy += -Yh11*(y12*w1xry+w1xry*y12);
				S1zz += -Yh11*(z12*w1xrz+w1xrz*z12);
				S1yz += -Yh11*(y12*w1xrz+w1xry*z12);
				S1zx += -Yh11*(z12*w1xrx+w1xrz*x12);
				S1xy += -Yh11*(x12*w1xry+w1xrx*y12);

				objects[n1_sph].f.x += Fx;			
				objects[n1_sph].f.y += Fy;
				objects[n1_sph].f.z += Fz;
				objects[n2_obj].f.x -= Fx;
				objects[n2_obj].f.y -= Fy;
				objects[n2_obj].f.z -= Fz;
				objects[n1_sph].t.x += T1x;
				objects[n1_sph].t.y += T1y;
				objects[n1_sph].t.z += T1z;
				objects[n1_sph].pf_lub.xx += S1xx;
				objects[n1_sph].pf_lub.yy += S1yy;
				objects[n1_sph].pf_lub.zz += S1zz;
				objects[n1_sph].pf_lub.yz += S1yz;
				objects[n1_sph].pf_lub.zx += S1zx;
				objects[n1_sph].pf_lub.xy += S1xy;


				psi = Xa/objects[n1_sph].mass;
				if (psi > Kappa)                                           /* List walls for implicit update */
				{
					objects[n1_sph].cl_list[0]++;
					objects[n1_sph].cl_list[objects[n1_sph].cl_list[0]] = n2_obj;
					if (objects[n1_sph].cl_list[0] > MAX_N ) fatal_err("too many neighbors in cluster", MAX_N);
				}
			}
		}
	}
}
