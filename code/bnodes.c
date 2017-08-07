/*  BNODES:  Node boundary conditions: Contains
		BNODES_SPH  Boundary conditions and forces for spheroids
		BNODES_WALL Boundary conditions and forces for walls  */
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
extern int    num_sph, num_obj;
extern int    num_proc, n_proc;


/*  BNODES_SPH_1: Computes forces and torques from boundary-node forces  */

void bnodes_sph_1 (struct object *objects, Float ***velcs_df, int **node_map)
{
	struct object *sph;
	double rx, ry, rz;
	double rbx, rby, rbz, fbx, fby, fbz;
	double fi, fci;
	int    x1, y1, z1, l1, x2, y2, z2, l2, xy2;
	int    node, link, n_sph, n_bnode, num_wall;

	num_wall = (num_obj - num_sph);                                     /* # of walls */

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{	
		sph  = &objects[n_sph];
		if (sph->num_bnode == 0)  continue;                             /* Skip if no boundary nodes */
		rx   = sph->r_map.x;
		ry   = sph->r_map.y;
		rz   = sph->r_map.z;
		if (range(box(rx, max_x), x_min, x_max))                        /* Particle center in processor domain */
		{
			sph->f.x -= sph->fluid*sph->c.x/sph->sum;
			sph->f.y -= sph->fluid*sph->c.y/sph->sum;
			sph->f.z -= sph->fluid*sph->c.z/sph->sum;
			sph->t.x -= sph->fluid*sph->rxc.x/sph->sum;
			sph->t.y -= sph->fluid*sph->rxc.y/sph->sum;
			sph->t.z -= sph->fluid*sph->rxc.z/sph->sum;
		}

		for (n_bnode = 0; n_bnode < sph->num_bnode; n_bnode++)
		{
			node = sph->nodes[n_bnode];                                 /* Boundary nodes */
			link = sph->links[n_bnode];                                 /* Boundary links */
			x1   = node/((max_y+2)*(max_z+2)) - 1;
			y1   = node%((max_y+2)*(max_z+2))/(max_z+2);
			z1   = node%(max_z+2);
			if (range (x1, 0, num_x))                                   /* Node in processor domain */
			{
				rbx  = x1 + 0.5 - rx + x_min;                           /* Coords (incl. box start) */
				rby  = y1 + 0.5 - ry;
				rbz  = z1 + 0.5 - rz;
				rbx  = n_image(rbx, max_x);                             /* Periodic BC's */
				rby  = n_image(rby, max_y);
				rbz  = n_image(rbz, max_z);
				fbx  = 0.0;
				fby  = 0.0;
				fbz  = 0.0;

				for (l1 = 0; l1 < Num_Dir; l1++)
				{
					if (link>>l1 & 1)                                   /* Boundary node */
					{
						x2  = x1 + c_x[l1];
						y2  = y1 + c_y[l1];
						z2  = z1 + c_z[l1];
						if ((num_wall > 0 && range(y2,0,max_y) == 0)    /* Shared node with wall */
								||  (num_wall > 2 && range(z2,0,max_z) == 0)
								||  (num_wall > 4 && range(x2+x_min,0,max_x) == 0))  continue;
						y2  = (y2+(max_y+1))%(max_y+1);
						z2  = (z2+(max_z+1))%(max_z+1);
						xy2 = (x2+1)*(max_y+2) + y2;

						if (node_map[xy2][z2]==0)          /* Fluid node */
						{
							l2   = l1 + (l1%2)*2 - 1;
							fi   = velcs_df[xy2][l2][z2];

							/* Completes particle-solid mass transfer across the boundary started in bnode_add */
							/* df = -2*a*ub dot c */
							fci  = -2*fi;
							fbx += fci*c_x[l1];                         /* Accumulate forces */
							fby += fci*c_y[l1];
							fbz += fci*c_z[l1];
						}
					}
				}
				sph->f.x += fbx;
				sph->f.y += fby;
				sph->f.z += fbz;
				sph->t.x += rby*fbz - rbz*fby;
				sph->t.y += rbz*fbx - rbx*fbz;
				sph->t.z += rbx*fby - rby*fbx;
			}
		}
	}
}


/*  BNODES_SPH_2: Updates LBE for boundary conditions  */

void bnodes_sph_2 (struct object *objects, Float ***velcs_df, int **node_map)
{
	struct object *sph;
	double rx, ry, rz, ux, uy, uz, wx, wy, wz;
	double rbx, rby, rbz, ubx, uby, ubz, fbx, fby, fbz;
	double rlx, rly, rlz;
	double fi, uci, fci, dfi;
	int    x1, y1, z1, l1, xy1, x2, y2, z2, l2, xy2;
	int    node, link, flag, n_sph, n_bnode, num_wall;

	num_wall = (num_obj - num_sph);                                     /* # of walls */

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{	
		sph = &objects[n_sph];
		dfi = sph->fluid/sph->sum;                                      /* Fluid from map update */
		sph->fluid = 0;
		if (sph->num_bnode == 0)  continue;                             /* Skip if no boundary nodes */
		rx  = sph->r_map.x;
		ry  = sph->r_map.y;
		rz  = sph->r_map.z;
		ux  = sph->u.x;
		uy  = sph->u.y;
		uz  = sph->u.z;
		wx  = sph->w.x;
		wy  = sph->w.y;
		wz  = sph->w.z;
		dfi-= (ux*sph->c.x + uy*sph->c.y + uz*sph->c.z + wx*sph->rxc.x + wy*sph->rxc.y + wz*sph->rxc.z)/sph->sum;

		for (n_bnode = 0; n_bnode < sph->num_bnode; n_bnode++)
		{
			node = sph->nodes[n_bnode];                                 /* Boundary nodes */
			link = sph->links[n_bnode];                                 /* Boundary links */
			x1   = node/((max_y+2)*(max_z+2)) - 1;
			y1   = node%((max_y+2)*(max_z+2))/(max_z+2);
			z1   = node%(max_z+2);
			xy1  = (x1+1)*(max_y+2) + y1;
			if (range (x1, 0, num_x))  flag = 1;                        /* Node in processor domain */
			else                       flag = 0;
			rbx  = x1 + 0.5 - rx + x_min;                               /* Coords (incl. box start) */
			rby  = y1 + 0.5 - ry;
			rbz  = z1 + 0.5 - rz;
			rbx  = n_image(rbx, max_x);                                 /* Periodic BC's */
			rby  = n_image(rby, max_y);
			rbz  = n_image(rbz, max_z);
			ubx  = ux + wy*rbz - wz*rby;
			uby  = uy + wz*rbx - wx*rbz;
			ubz  = uz + wx*rby - wy*rbx;

			for (l1 = 0; l1 < Num_Dir; l1++)
			{
				if (link>>l1 & 1)                                       /* Boundary node */
				{
					x2  = x1 + c_x[l1];
					y2  = y1 + c_y[l1];
					z2  = z1 + c_z[l1];
					if ((num_wall > 0 && range(y2,0,max_y) == 0)        /* Shared node with wall */
							||  (num_wall > 2 && range(z2,0,max_z) == 0)
							||  (num_wall > 4 && range(x2+x_min,0,max_x) == 0))  continue;
					y2  = (y2+(max_y+2))%(max_y+2);
					z2  = (z2+(max_z+2))%(max_z+2);
					xy2 = (x2+1)*(max_y+2) + y2;

					if (node_map[xy2][z2]==0)              /* Fluid node */
					{
						l2   = l1 + (l1%2)*2 - 1;
						fi   = velcs_df[xy2][l2][z2];
						uci  = 2*fac[l1]*(ubx*c_x[l1] + uby*c_y[l1] + ubz*c_z[l1])/CS2;
						uci += 2*fac[l1]*dfi/CS2;                       /* Correct for mass transfer across fluid-solid boundary */
						fci  = -(2*fi + uci)*flag;
						fbx  = fci*c_x[l1];                             /* Accumulate forces */
						fby  = fci*c_y[l1];
						fbz  = fci*c_z[l1];
						rlx = rbx + 0.5*c_x[l1];
						rly = rby + 0.5*c_y[l1];
						rlz = rbz + 0.5*c_z[l1];

						/* update particle-fluid stresses */
						sph->pf_str.xx +=  rlx*fbx;
						sph->pf_str.yy +=  rly*fby;
						sph->pf_str.zz +=  rlz*fbz;
						sph->pf_str.yz += (rly*fbz + rlz*fby)/2.0;
						sph->pf_str.zx += (rlz*fbx + rlx*fbz)/2.0;
						sph->pf_str.xy += (rlx*fby + rly*fbx)/2.0;
						velcs_df[xy1][l1][z1] = fi + uci;
					}
				}
			}
		}
	}
}


/*  BNODES_WALL: Updates LBE for boundary conditions at plane walls  */

void bnodes_wall (struct object *objects, Float ***velcs_df, int **node_map)
{
	double ux, uy, uz, fx, fy, fz;
	double fi, uci, fci;
	int    i1, j1, max_i, max_j;
	int    x1, y1, z1, l1, xy1, x2, y2, z2, l2, xy2;
	int    n, n_obj, n_wall, flag;

	for (n_obj = num_sph; n_obj < num_obj; n_obj++)
		for (n = 0;  n < Num_Dir_X; n++)
		{
			n_wall = n_obj - num_sph;
			ux = objects[n_obj].u.x;
			uy = objects[n_obj].u.y;
			uz = objects[n_obj].u.z;
			fx = fy = fz = 0.0;

			switch (n_wall)
			{
				case 0:                                                     /* Y- wall */
					max_i = max_z + 1;
					max_j = num_x + 1;
					l1  = y_m[n];
					l2  = y_p[n];
					break;
				case 1:                                                     /* Y+ wall */
					max_i = max_z + 1;
					max_j = num_x + 1;
					l1  = y_p[n];
					l2  = y_m[n];
					break;
				case 2:                                                     /* Z- wall */
					max_i = num_x + 1;
					max_j = max_y + 1;
					l1  = z_m[n];
					l2  = z_p[n];
					break;
				case 3:                                                     /* Z+ wall */
					max_i = num_x + 1;
					max_j = max_y + 1;
					l1  = z_p[n];
					l2  = z_m[n];
					break;
				case 4:                                                     /* X- wall */
					if (n_proc != 0)  continue;                             /* First proc only */
					max_i = max_y + 1;
					max_j = max_z + 1;
					l1  = x_m[n];
					l2  = x_p[n];
					break;
				case 5:                                                     /* X+ wall */
					if (n_proc != num_proc-1)  continue;                    /* Last proc only */
					max_i = max_y + 1;
					max_j = max_z + 1;
					l1  = x_p[n];
					l2  = x_m[n];
					break;
			}
			for(i1 = 1; i1 <  max_i; i1++)
				for(j1 = 1; j1 <  max_j; j1++)
				{
					switch (n_wall)
					{
						case 0:                                                 /* Y- wall */
							x1 = j1;
							y1 = 1;
							z1 = i1;
							break;
						case 1:                                                 /* Y+ wall */
							x1 = j1;
							y1 = max_y;
							z1 = i1;
							break;
						case 2:                                                 /* Z- wall */
							x1 = i1;
							y1 = j1;
							z1 = 1;
							break;
						case 3:                                                 /* Z+ wall */
							x1 = i1;
							y1 = j1;
							z1 = max_z;
							break;
						case 4:                                                 /* X- wall; first proc only */
							x1 = 1;
							y1 = i1;
							z1 = j1;
							break;
						case 5:                                                 /* X+ wall; last proc only */
							x1 = num_x;
							y1 = i1;
							z1 = j1;
							break;
					}
					x2  = x1 + c_x[l1];
					y2  = y1 + c_y[l1];
					z2  = z1 + c_z[l1];
					xy1 = x1*(max_y+2) + y1;
					xy2 = x2*(max_y+2) + y2;

					flag = range(x1, 1, num_x+1);       /* Node in proc domain (flag = 1) */
					if (flag || range(x2, 1, num_x+1)){  /* Other node in proc domain */
						if (n_wall > 1) {
							if ((y1 == 0)       && (c_y[l1] == -1))  continue;  /* Skip XY & YZ corners */
							if ((y1 == max_y) && (c_y[l1] ==  1))  continue;
						}
						if (n_wall > 3) {
							if ((z1 == 0)       && (c_z[l1] == -1))  continue;  /* Skip XZ corners */
							if ((z1 == max_z) && (c_z[l1] ==  1))  continue;
						}

						/* momentum transfer across solid-fluid boundary */
						/* uci = u_b dot c */
						uci = (ux*c_x[l1] + uy*c_y[l1] + uz*c_z[l1])*fac[l1]*2/CS2;
						if((num_obj-num_sph == 4) && (n_wall <= 1)) {            /* 4-wall elongational flow */
							if((z1 == 1)     && (c_z[l1] == -1))  uci = 0;   /* Set corner velocity to zero */
							if((z1 == max_z) && (c_z[l1] ==  1))  uci = 0;
						}

						if(node_map[xy1][z1]==0)        /* Fluid node */
							fi  = velcs_df[xy1][l1][z1];
						else {                          /* Shared node */ 
							fi  = 0.0;
							uci = 0.0;
						}
						velcs_df[xy2][l2][z2] = fi - uci;
					}
					fx += fci*c_x[l1];
					fy += fci*c_y[l1];
					fz += fci*c_z[l1];
				}

			objects[n_obj].f.x += fx;
			objects[n_obj].f.y += fy;
			objects[n_obj].f.z += fz;
		}
}
