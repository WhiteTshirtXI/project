/*  N_LIST: Construct near-neighbor list  */
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

void n_list (struct object *objects)
{
	extern int  max_x, max_y, max_z;
	extern int  num_sph;
	static int  cell_list[MAX_L][MAX_N+1], flag=0;
	double x12, y12, z12, r12, r_max, sig, cut;
	int nlx, nly, nlz;
	int lx1, ly1, lz1, l1, n1, lx2, ly2, lz2, l2, n2;
	int n1_sph, n2_sph;


	for (n1_sph = 0, r_max = 0.0; n1_sph < num_sph; n1_sph++)          /* Largest radius */
	{
		r_max = max(r_max, objects[n1_sph].r_a);
		objects[n1_sph].r_lst.x = objects[n1_sph].r.x;                 /* Store positions */
		objects[n1_sph].r_lst.y = objects[n1_sph].r.y;
		objects[n1_sph].r_lst.z = objects[n1_sph].r.z;
		objects[n1_sph].list[0] = 0;
	}

	cut = 2.0*r_max + 2.0;                                               /* Set cell sizes */
	nlx = max_x/cut; nly = max_y/cut; nlz = max_z/cut;	               /* Set cell numbers */


	/*  Calculate neighbor lists  */

	if (min(nlx, min(nly, nlz)) < 3)                                   /* Too few cells: use direct sum */
	{
		for (n1_sph = 1; n1_sph < num_sph; n1_sph++)
			for (n2_sph = 0; n2_sph < n1_sph;  n2_sph++)
			{
				x12 = n_image(objects[n1_sph].r.x - objects[n2_sph].r.x, max_x);
				y12 = n_image(objects[n1_sph].r.y - objects[n2_sph].r.y, max_y);
				z12 = n_image(objects[n1_sph].r.z - objects[n2_sph].r.z, max_z);
				r12 = x12*x12 + y12*y12 + z12*z12;
				sig = objects[n1_sph].r_a + objects[n2_sph].r_a;
				cut = sig + 2.0;
				if (r12 <= cut*cut)
				{
					if (r12 - sig*sig < -Tol)  {
						warning ("overlap in n_list");
						printf("Particle %d at (%g %g %g) and Particle %d at (%g %g %g) overlap, r12=%le sig2=%le\n ", n1_sph, objects[n1_sph].r.x, objects[n1_sph].r.y, objects[n1_sph].r.z, n2_sph, objects[n2_sph].r.x, objects[n2_sph].r.y, objects[n2_sph].r.z, r12, sig*sig);
					}
					objects[n1_sph].list[0]++;
					objects[n2_sph].list[0]++;
					if (objects[n1_sph].list[0] > MAX_N)
						fatal_err("too many neighbors", MAX_N);
					if (objects[n2_sph].list[0] > MAX_N)
						fatal_err("too many neighbors", MAX_N);
					objects[n1_sph].list[objects[n1_sph].list[0]] = n2_sph;
					objects[n2_sph].list[objects[n2_sph].list[0]] = n1_sph;
				}
			}
	}
	else                                                               /* Use cell structure */
	{

		while (nlx*nly*nlz > MAX_L)
		{
			if (flag == 0)  warning ("too many cells; doubling cell size");
			nlx /= 2;  nly /= 2;  nlz /=2;
		}
		flag = 1;                                                      /* Turn warning off */

		for (lx1 = 0; lx1 < nlx; lx1++)
			for (ly1 = 0; ly1 < nly; ly1++)
				for (lz1 = 0; lz1 < nlz; lz1++)
				{
					l1  = lx1*nly*nlz + ly1*nlz + lz1;
					cell_list[l1][0] = 0;
				}

		for (n1_sph = 0; n1_sph < num_sph; n1_sph++)
		{
			lx1 = (int) (box(objects[n1_sph].r.x/max_x, 1)*nlx);
			ly1 = (int) (box(objects[n1_sph].r.y/max_y, 1)*nly);
			lz1 = (int) (box(objects[n1_sph].r.z/max_z, 1)*nlz);
			l1  = lx1*nly*nlz + ly1*nlz + lz1;
			cell_list[l1][0]++;
			if (cell_list[l1][0] > MAX_N)
				fatal_err("too many cell entries:", -1);
			cell_list[l1][cell_list[l1][0]] = n1_sph;
			objects[n1_sph].list[0] = 0;
		}

		for (lx1 = 0; lx1 < nlx; lx1++)
			for (ly1 = 0; ly1 < nly; ly1++)
				for (lz1 = 0; lz1 < nlz; lz1++)
					for (lx2 = lx1-1; lx2 <= lx1+1; lx2++)
						for (ly2 = ly1-1; ly2 <= ly1+1; ly2++)
							for (lz2 = lz1-1; lz2 <= lz1+1; lz2++)
							{
								l1  = lx1*nly*nlz + ly1*nlz + lz1;
								l2  = ((lx2+nlx)%nlx)*nly*nlz + ((ly2+nly)%nly)*nlz + (lz2+nlz)%nlz;
								for (n1 = 1; n1 <= cell_list[l1][0]; n1++)
									for (n2 = 1; n2 <= cell_list[l2][0]; n2++)
									{
										n1_sph  = cell_list[l1][n1];
										n2_sph  = cell_list[l2][n2];
										if (n2_sph < n1_sph)
										{
											x12 = n_image(objects[n1_sph].r.x - objects[n2_sph].r.x, max_x);
											y12 = n_image(objects[n1_sph].r.y - objects[n2_sph].r.y, max_y);
											z12 = n_image(objects[n1_sph].r.z - objects[n2_sph].r.z, max_z);
											r12 = x12*x12 + y12*y12 + z12*z12;
											sig = objects[n1_sph].r_a + objects[n2_sph].r_a;
											cut = sig + 2.0;
											if (r12 <= cut*cut)
											{
												if (r12 - sig*sig < -Tol)  {
													warning ("overlap in n_list");
													printf("Particle %d at (%g %g %g) and Particle %d at (%g %g %g) overlap, r12=%le cut2=%le\n ", n1_sph, objects[n1_sph].r.x, objects[n1_sph].r.y, objects[n1_sph].r.z, n2_sph, objects[n2_sph].r.x, objects[n2_sph].r.y, objects[n2_sph].r.z, r12, cut*cut);
												}
												objects[n1_sph].list[0]++;
												objects[n2_sph].list[0]++;
												if (objects[n1_sph].list[0] > MAX_N)  fatal_err("too many neighbors", MAX_N);
												if (objects[n2_sph].list[0] > MAX_N)  fatal_err("too many neighbors", MAX_N);
												objects[n1_sph].list[objects[n1_sph].list[0]] = n2_sph;
												objects[n2_sph].list[objects[n2_sph].list[0]] = n1_sph;
											}
										}
									}
							}
	}
}


void n_list_mon(struct sphere_param *sphere_pm, struct monomer *mon)
{
	extern int  max_x, max_y, max_z, wall_flag;
	static int  cell_list[MAX_L][MAX_N+1], flag=0;
	double x12, y12, z12, r12, r_max, sig, cut;
	int nlx, nly, nlz;
	int lx1, ly1, lz1, l1, n1, lx2, ly2, lz2, l2, n2;
	int n1_sph, n2_sph;
	int num_beads = sphere_pm->num_beads;

	for(n1_sph = 0, r_max = 0.0; n1_sph < num_beads; n1_sph++) {  // Largest radius
		r_max = max(r_max, mon[n1_sph].radius);
		mon[n1_sph].pos_lst[0] = mon[n1_sph].pos[0];      // Store positions
		mon[n1_sph].pos_lst[1] = mon[n1_sph].pos[1];
		mon[n1_sph].pos_lst[2] = mon[n1_sph].pos[2];
		mon[n1_sph].list[0] = 0;
	}
	//cut = sphere_pm->evcutoff*r_max+2.0;              // Set cell sizes
  cut=4.; // Modification 20160317
 	nlx = max_x/cut; nly = max_y/cut; nlz = max_z/cut;  // Set cell numbers

	//  Calculate neighbor lists
	if(min(nlx, min(nly, nlz)) < 3)  // Too few cells: use direct sum 
	{
		for(n1_sph = 1; n1_sph < num_beads; n1_sph++) {
			for(n2_sph = 0; n2_sph < n1_sph;  n2_sph++)
			{
				x12 = n_image(mon[n1_sph].pos[0] - mon[n2_sph].pos[0], max_x);
				y12 = n_image(mon[n1_sph].pos[1] - mon[n2_sph].pos[1], max_y);
				z12 = n_image(mon[n1_sph].pos[2] - mon[n2_sph].pos[2], max_z);
				r12 = x12*x12 + y12*y12 + z12*z12;
				sig = (mon[n1_sph].radius + mon[n2_sph].radius)/2.0;
				cut = sig*sphere_pm->evcutoff;
				if(r12 <= cut*cut)
				{
					mon[n1_sph].list[0]++;
					mon[n2_sph].list[0]++;
					if(mon[n1_sph].list[0] > MAX_N)
					  fatal_err("too many neighbors", MAX_N);
					if(mon[n2_sph].list[0] > MAX_N)
					  fatal_err("too many neighbors", MAX_N);
					mon[n1_sph].list[mon[n1_sph].list[0]] = n2_sph;
					mon[n2_sph].list[mon[n2_sph].list[0]] = n1_sph;
				}
				/*	printf("monomer %d neighbor %d = %d\n", n1_sph, mon[n1_sph].list[0], n2_sph); */
			}
		}
	}
	else  // Use cell structure 
	{
		while(nlx*nly*nlz > MAX_L) {
			if(flag == 0)  
        warning ("too many cells; doubling cell size");
			nlx /= 2;    nly /= 2;    nlz /= 2;
		}
		flag = 1;  // Turn warning off

		for(lx1 = 0; lx1 < nlx; lx1++)
			for(ly1 = 0; ly1 < nly; ly1++)
				for(lz1 = 0; lz1 < nlz; lz1++) {
					l1  = lx1*nly*nlz + ly1*nlz + lz1;
					cell_list[l1][0] = 0;
				}

		for(n1_sph = 0; n1_sph < num_beads; n1_sph++)
		{
			lx1 = (int) (box(mon[n1_sph].pos[0]/max_x, 1)*nlx);
			ly1 = (int) (box(mon[n1_sph].pos[1]/max_y, 1)*nly);
			lz1 = (int) (box(mon[n1_sph].pos[2]/max_z, 1)*nlz);
			l1  = lx1*nly*nlz + ly1*nlz + lz1;
			cell_list[l1][0]++;
			if (cell_list[l1][0] > MAX_N)
				fatal_err("too many cell entries:", -1);
			cell_list[l1][cell_list[l1][0]] = n1_sph;
			mon[n1_sph].list[0] = 0;
		}

		for(lx1 = 0; lx1 < nlx; lx1++)
			for(ly1 = 0; ly1 < nly; ly1++)
				for(lz1 = 0; lz1 < nlz; lz1++)
					for(lx2 = lx1-1; lx2 <= lx1+1; lx2++)
						for(ly2 = ly1-1; ly2 <= ly1+1; ly2++)
							for(lz2 = lz1-1; lz2 <= lz1+1; lz2++)
							{
								l1  = lx1*nly*nlz + ly1*nlz + lz1;
								l2  = ((lx2+nlx)%nlx)*nly*nlz + ((ly2+nly)%nly)*nlz + (lz2+nlz)%nlz;
								for(n1 = 1; n1 <= cell_list[l1][0]; n1++)
									for(n2 = 1; n2 <= cell_list[l2][0]; n2++)
									{
										n1_sph  = cell_list[l1][n1];
										n2_sph  = cell_list[l2][n2];
										if(n2_sph < n1_sph)
										{
											x12 = n_image(mon[n1_sph].pos[0] - mon[n2_sph].pos[0], max_x);
											y12 = n_image(mon[n1_sph].pos[1] - mon[n2_sph].pos[1], max_y);
											z12 = n_image(mon[n1_sph].pos[2] - mon[n2_sph].pos[2], max_z);
											if(wall_flag > 0)
												y12 = mon[n1_sph].pos[1]-mon[n2_sph].pos[1];
											if(wall_flag > 1)
												z12 = mon[n1_sph].pos[2]-mon[n2_sph].pos[2];
											r12 = x12*x12 + y12*y12 + z12*z12;
											sig = (mon[n1_sph].radius+mon[n2_sph].radius)/2.0;
											//cut = 2.3;
                      //cut=3.8;
                      cut=3.;
//                      cut = (sphere_pm->cutoff_LJ+0.5) * sphere_pm->range_LJ; 
                      //cut = 3.;  //20170802
											//if(r12 <= cut*cut)  // Modification 20170318!!!
                      if(sqrt(r12) <= cut) 
                      {
												mon[n1_sph].list[0]++;
												mon[n2_sph].list[0]++;
												if(mon[n1_sph].list[0] > MAX_N)
													fatal_err("too many neighbors", MAX_N);
												if(mon[n2_sph].list[0] > MAX_N)
													fatal_err("too many neighbors", MAX_N);
												mon[n1_sph].list[mon[n1_sph].list[0]] = n2_sph;
												mon[n2_sph].list[mon[n2_sph].list[0]] = n1_sph;
											}
										}
									}
							}
	}
}

