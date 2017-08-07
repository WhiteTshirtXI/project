/* CLUSTER_SINGULAR_MATRIX: calculate implicit velocity for clusters */
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

void cluster_update (struct object *objects, struct cluster *clusters, double dt, int n_cl) 
{
    double 	   **matrix, *tmp_p, *u;
	int		   n, n_sph, matrix_size;

	matrix_size = clusters[n_cl].list[0];
	if (matrix_size > MAX_Cl)  fatal_err ("Cluster matrix too large", matrix_size);

	matrix_size *= 3;
    matrix = (double **) calloc(matrix_size, sizeof(*matrix));
    tmp_p  = (double *)  calloc(matrix_size*matrix_size, sizeof(*tmp_p));
	u = (double *)       calloc(matrix_size, sizeof(*u));
	if (matrix == 0)  fatal_err ("cannot allocate matrix", matrix_size);
	if (tmp_p == 0)   fatal_err ("cannot allocate tmp_p", matrix_size);
	if (u == 0)	      fatal_err ("cannot allocate u", matrix_size);
    for (n = 0; n < matrix_size; n++)
    {
        matrix[n] = tmp_p;
        tmp_p    += matrix_size;
    }

	for (n = 0; n < clusters[n_cl].list[0]; n++)
	{
		n_sph = clusters[n_cl].list[n+1]; 
        u[n*3  ] = objects[n_sph].u.x;
        u[n*3+1] = objects[n_sph].u.y;
        u[n*3+2] = objects[n_sph].u.z;
	}

   	cluster_matrix (objects, clusters, matrix, n_cl, dt);
    cj_grad (matrix, u, matrix_size);

	for (n = 0; n < clusters[n_cl].list[0]; n++)
	{
		n_sph = clusters[n_cl].list[n+1]; 
		objects[n_sph].u.x  = u[n*3  ];
		objects[n_sph].u.y  = u[n*3+1];
		objects[n_sph].u.z  = u[n*3+2];
	}

	free (u);
	free (matrix[0]);
	free (matrix);
 }

void cluster_matrix (struct object *objects, struct cluster *clusters, double **matrix, int n_cl, double dt)
{
    extern double tau;
    extern double lub_N;
	extern int    max_x, max_y, max_z, num_obj;
    extern int    num_sph, num_obj;
    double  nu, dr, rad, cut;
    double  alpha, beta1, beta4, Xa, Xa_max;
    double  x12, y12, z12, r12;
	int     *index;
    int     n1_sph, n2_obj, n, n1, n2, matrix_size;

    nu = (tau - 0.5)/3.0;

	index = (int *)  calloc(num_sph, sizeof(int));                 /* Create temporary indexing array */
	if (index == 0)  fatal_err("cannot allocate index", num_sph);
	for (n = 0; n < clusters[n_cl].list[0]; n++)                     
		index[clusters[n_cl].list[n+1]] = n;

    for (n1 = 0; n1 < clusters[n_cl].list[0]; n1++)
    {
		n1_sph = clusters[n_cl].list[n1+1];

        for (n = 0; n < objects[n1_sph].cl_list[0]; n++)
        {
       		n2_obj = objects[n1_sph].cl_list[n+1];

			if (n2_obj < num_sph)
			{
			    x12    = n_image(objects[n1_sph].r.x - objects[n2_obj].r.x, max_x);
   		     	y12    = n_image(objects[n1_sph].r.y - objects[n2_obj].r.y, max_y);
       		    z12    = n_image(objects[n1_sph].r.z - objects[n2_obj].r.z, max_z);
				r12    = sqrt(x12*x12 + y12*y12 + z12*z12);
       	    	rad    = (objects[n1_sph].r_a + objects[n2_obj].r_a)/2.0;
       	    	dr     = r12 - 2.0*rad;
			}
			else 
			{
       	    	x12 = (objects[n1_sph].r.x - objects[n2_obj].r.x)*objects[n2_obj].e.x*objects[n2_obj].e.x;
       	    	y12 = (objects[n1_sph].r.y - objects[n2_obj].r.y)*objects[n2_obj].e.y*objects[n2_obj].e.y;
      	     	z12 = (objects[n1_sph].r.z - objects[n2_obj].r.z)*objects[n2_obj].e.z*objects[n2_obj].e.z;
       	    	r12 = sqrt(x12*x12 + y12*y12 + z12*z12);
      	     	rad = objects[n1_sph].r_a;
				dr  = r12 - rad;
			}

        	alpha = dt*Pi*Rho_Fl*nu*rad;
       		cut = (objects[n1_sph].lub_cut + objects[n2_obj].lub_cut)/2.0;
			dr = max(dr, cut);
						
			if (n2_obj < num_sph)
			{
				beta1 = objects[n2_obj].r_a/objects[n1_sph].r_a;
				beta4 = pow((1 + beta1), 4.0);
				Xa    =  4.0*beta1*beta1/beta4;
			}
			else
				Xa  = 1.0;

			if (objects[n2_obj].mass_flag)  Xa_max = Kappa*min(objects[n1_sph].mass, objects[n2_obj].mass);
			else                            Xa_max = Kappa*objects[n1_sph].mass;

			Xa *= 6.0*alpha*rad/dr;
			Xa -= Xa_max;
			Xa /= objects[n1_sph].mass;

			x12   /= r12;
			y12   /= r12;
        	z12   /= r12;

			if (objects[n2_obj].mass_flag)
			{
           		n2 = index[n2_obj];
            	matrix[n1*3  ][n1*3  ] +=  Xa*x12*x12;
            	matrix[n1*3  ][n1*3+1] +=  Xa*x12*y12;
            	matrix[n1*3  ][n1*3+2] +=  Xa*x12*z12;
            	matrix[n1*3  ][n2*3  ]  = -Xa*x12*x12;
            	matrix[n1*3  ][n2*3+1]  = -Xa*x12*y12;
            	matrix[n1*3  ][n2*3+2]  = -Xa*x12*z12;
            	matrix[n1*3+1][n1*3  ] +=  Xa*y12*x12;
            	matrix[n1*3+1][n1*3+1] +=  Xa*y12*y12;
            	matrix[n1*3+1][n1*3+2] +=  Xa*y12*z12;
            	matrix[n1*3+1][n2*3  ]  = -Xa*y12*x12;
            	matrix[n1*3+1][n2*3+1]  = -Xa*y12*y12;
            	matrix[n1*3+1][n2*3+2]  = -Xa*y12*z12;
            	matrix[n1*3+2][n1*3  ] +=  Xa*z12*x12;
            	matrix[n1*3+2][n1*3+1] +=  Xa*z12*y12;
            	matrix[n1*3+2][n1*3+2] +=  Xa*z12*z12;
            	matrix[n1*3+2][n2*3  ]  = -Xa*z12*x12;
            	matrix[n1*3+2][n2*3+1]  = -Xa*z12*y12;
            	matrix[n1*3+2][n2*3+2]  = -Xa*z12*z12;
			}
			else
			{
        	    matrix[n1*3  ][n1*3  ] +=  Xa*x12*x12;
				matrix[n1*3  ][n1*3+1] +=  Xa*x12*y12;
				matrix[n1*3  ][n1*3+2] +=  Xa*x12*z12;
				matrix[n1*3+1][n1*3  ] +=  Xa*y12*x12;
				matrix[n1*3+1][n1*3+1] +=  Xa*y12*y12;
				matrix[n1*3+1][n1*3+2] +=  Xa*y12*z12;
				matrix[n1*3+2][n1*3  ] +=  Xa*z12*x12;
				matrix[n1*3+2][n1*3+1] +=  Xa*z12*y12;
				matrix[n1*3+2][n1*3+2] +=  Xa*z12*z12;
			}
		}
    }
 
    matrix_size = 3*clusters[n_cl].list[0];
    for (n = 0; n < matrix_size; n++)                                         /* Add diagonal element */
        matrix[n][n] ++;
 
	free (index);
}
