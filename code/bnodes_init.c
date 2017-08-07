/*  BNODES_INIT:  Initializes node boundary conditions: Contains
		BNODES_INIT  Initialize boundary nodes for spheroids
		BNODES_ADD   Creates fluid in old boundary node regions
		BNODES_DEL   Deletes fluid in new boundary node regions
		BNODES_MOM   Calculate boundary nodes moments for velocity update
		MATRIX_INV   Calculate mobility matrices
		GAUSS-JORDAN Gauss-Jordan inversion  */
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


/*  BNODES_INIT: Create b_node map in processor domain  */

void bnodes_init (struct object *objects, int **node_map)
{
	struct object *sph;
	double rx, ry, rz;
	double del_x, del_y, del_z;
	int    size_xy, size_z, x, y, z, n;
	int    x1, y1, z1, l1, xy1, xp, bnode_flag;
	int    link, flag_int, flag_ext, n_sph, num_wall;

	num_wall = (num_obj - num_sph);                                         /* # of walls*/

	size_xy  = (num_x+2)*(max_y+2);
	size_z   = max_z+2;

	for (n = 0; n < size_xy; n++)     /* Initialize node_map */
		for (z = 0; z < size_z; z++)
			node_map[n][z] = 0;

	n=size_xy/2;

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{
		sph = &objects[n_sph];
		rx  = box(sph->r.x, max_x);   /* New boxed coords */
		ry  = box(sph->r.y, max_y);
		rz  = box(sph->r.z, max_z);
		sph->r_map.x = sph->r.x;      /* Origin of new boundary node map */
		sph->r_map.y = sph->r.y;
		sph->r_map.z = sph->r.z;
		sph->vol = 0;                 /* Zero interior fluid volume */

		if (sph->nodes != 0)          /* free unused space */
		{
			free (sph->nodes);
			free (sph->links);
			sph->nodes = 0;
			sph->num_bnode = 0;
		}

		for (xp = -max_x, bnode_flag=0; xp <= max_x; xp += max_x)           /* Check processor domain */
			if (range(rx+xp, x_min-sph->size-1,x_max+sph->size+1))  bnode_flag = 1;

		if (bnode_flag)
		{
			sph->nodes = (int *) calloc(sph->max_bnode, sizeof(int));       /* Allocate space to nodes */
			sph->links = (int *) calloc(sph->max_bnode, sizeof(int));       /* Allocate space to links */
			if (sph->nodes == 0 || sph->links == 0) fatal_err("cannot allocate boundary nodes", -1);
			for (x = -sph->size; x <= sph->size; x++)
				for (y = -sph->size; y <= sph->size; y++)
					for (z = -sph->size; z <= sph->size; z++)
					{
						x1 = x + (int)rx;
						y1 = y + (int)ry;
						z1 = z + (int)rz;
						del_x  = x1 + 0.5 - rx;                                     /* Node coords = n+1/2 */ 
						del_y  = y1 + 0.5 - ry;
						del_z  = z1 + 0.5 - rz;
						del_x  = n_image(del_x, max_x);                             /* Periodic BC's */
						del_y  = n_image(del_y, max_y);
						del_z  = n_image(del_z, max_z);
						if (sph->func_ptr(sph, del_x, del_y, del_z)) continue;      /* Node outside particle */

						for (xp = x1-max_x; xp <= x1+max_x; xp += max_x)            /* Loop over periodic images in x */
						{
							if (range(xp, x_min-1, x_max+1) == 0)       continue;   /* Include nodes in x borders */
							if (num_wall > 4 && (xp==-1 || xp==max_x))  continue;   /* X-walls: no exterior nodes */
							x1  =  xp - x_min;                                      /* Local coordinate node numbers*/
							y1  = (y1+(max_y+1))%(max_y+1);
							z1  = (z1+(max_z+1))%(max_z+1);
							xy1 = (x1+1)*(max_y+2) + y1;
							flag_int = range(xp, x_min, x_max);                     /* Interior node in box (flag = 1) */
							sph->vol += flag_int;                                   /* Increment particle volume */
							node_map[xy1][z1] = 1;                /* Set solid bit in node map */

							for (l1 = 0, link = 0; l1 < Num_Dir; l1++)
							{
								flag_ext = range(xp+c_x[l1], x_min, x_max);         /* Exterior node in box (flag = 1) */
								if ((flag_int || flag_ext) == 0)  continue;         /* Include exterior nodes */
								flag_ext = sph->func_ptr(sph, del_x+c_x[l1], del_y+c_y[l1], del_z+c_z[l1]);
								if (flag_ext)  link |= 1 << l1;                     /* Set link bit */
							}
							if (link != 0)
							{
								sph->nodes[sph->num_bnode] = xy1*(max_z+2) + z1;
								sph->links[sph->num_bnode] = link;
								sph->num_bnode++;
								if (sph->num_bnode > sph->max_bnode)  fatal_err("too many b-nodes,", sph->num_bnode);
							}
						}
					}
		}
	}	
}


/*  BNODES_MOM: Boundary node moments for implicit velocity update  */
/*  Calculate the friction matrix of the particle */

void bnodes_mom (struct object *objects, int **node_map)
{
	struct object  *sph;
	double rx, ry, rz, a;
	double cx, cy, cz, rbx, rby, rbz, rxcx, rxcy, rxcz;
	int    x1, y1, z1, l1, x2, y2, z2, xy2;
	int    node, link, n_sph, n_bnode, num_wall;

	num_wall = (num_obj - num_sph);                                         /* # of walls*/

	/* add contributions to moment due to mass transfer across fluid-solid boundary */

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{	
		sph = &objects[n_sph];
		rx  = box(sph->r_map.x, max_x);
		ry  = box(sph->r_map.y, max_y);
		rz  = box(sph->r_map.z, max_z);
		sph->sum    = 0;                                                    /* Zero moments */
		sph->c.x    = 0;  sph->c.y    = 0;  sph->c.z    = 0;
		sph->rxc.x  = 0;  sph->rxc.y  = 0;  sph->rxc.z  = 0;
		sph->ztt.xx = 0;  sph->ztt.yy = 0;  sph->ztt.zz = 0;
		sph->ztt.yz = 0;  sph->ztt.zx = 0;  sph->ztt.xy = 0;
		sph->zrr.xx = 0;  sph->zrr.yy = 0;  sph->zrr.zz = 0;
		sph->zrr.yz = 0;  sph->zrr.zx = 0;  sph->zrr.xy = 0;
		sph->ztr.xx = 0;  sph->ztr.yy = 0;  sph->ztr.zz = 0;
		sph->ztr.yz = 0;  sph->ztr.zx = 0;  sph->ztr.xy = 0;
		sph->ztr.zy = 0;  sph->ztr.xz = 0;  sph->ztr.yx = 0;

		for (n_bnode = 0; n_bnode < sph->num_bnode; n_bnode++)
		{
			node = sph->nodes[n_bnode];                                     /* Boundary inside nodes */
			link = sph->links[n_bnode];                                     /* Boundary links */
			x1   = node/((max_y+2)*(max_z+2)) - 1;                                  /* Get coordinates from node index */
			y1   = node%((max_y+2)*(max_z+2))/(max_z+2);
			z1   = node%(max_z+2);
			if (range (x1, 0, num_x))                                       /* Node in processor domain */
			{
				rbx  = x1 + 0.5 - rx + x_min;                               /* Coords (incl. box start) */
				rby  = y1 + 0.5 - ry;
				rbz  = z1 + 0.5 - rz;
				rbx  = n_image(rbx, max_x);                                 /* Periodic BC's */
				rby  = n_image(rby, max_y);
				rbz  = n_image(rbz, max_z);

				for (l1 = 0; l1 < Num_Dir; l1++)
				{
					if (link>>l1 & 1)                                       /* Boundary node */
					{
						x2  = x1 + c_x[l1];
						y2  = y1 + c_y[l1];
						z2  = z1 + c_z[l1];
						if ((num_wall > 0 && range(y2,0,max_y+1) == 0)        /* Shared node with wall */
								||  (num_wall > 2 && range(z2,0,max_z+1) == 0)
								||  (num_wall > 4 && range(x2+x_min,0,max_x) == 0))  continue;
						y2  = (y2+(max_y+1))%(max_y+1);
						z2  = (z2+(max_z+1))%(max_z+1);
						xy2 = (x2+1)*(max_y+2) + y2;

						if (node_map[xy2][z2]==0)              /* Fluid node */
						{
							cx   = c_x[l1];
							cy   = c_y[l1];
							cz   = c_z[l1];
							rxcx = rby*cz - rbz*cy;                         /* r X c */
							rxcy = rbz*cx - rbx*cz;
							rxcz = rbx*cy - rby*cx;
							a    = 2*fac[l1]/CS2;

							/* boundary node contribution to the friction matrix */
							sph->sum    += a;
							sph->c.x    += a*cx;         sph->c.y    += a*cy;         sph->c.z    += a*cz;
							sph->rxc.x  += a*rxcx;       sph->rxc.y  += a*rxcy;       sph->rxc.z  += a*rxcz;
							sph->ztt.xx += a*cx*cx;      sph->ztt.yy += a*cy*cy;      sph->ztt.zz += a*cz*cz;
							sph->ztt.yz += a*cy*cz;      sph->ztt.zx += a*cz*cx;      sph->ztt.xy += a*cx*cy;
							sph->zrr.xx += a*rxcx*rxcx;  sph->zrr.yy += a*rxcy*rxcy;  sph->zrr.zz += a*rxcz*rxcz;
							sph->zrr.yz += a*rxcy*rxcz;  sph->zrr.zx += a*rxcz*rxcx;  sph->zrr.xy += a*rxcx*rxcy;
							sph->ztr.xx += a*cx*rxcx;    sph->ztr.yy += a*cy*rxcy;    sph->ztr.zz += a*cz*rxcz;
							sph->ztr.yz += a*cy*rxcz;    sph->ztr.zx += a*cz*rxcx;    sph->ztr.xy += a*cx*rxcy;
							sph->ztr.zy += a*cz*rxcy;    sph->ztr.xz += a*cx*rxcz;    sph->ztr.yx += a*cy*rxcx;
						}
					}
				}
			}
		}
	}

	if(num_sph > 0)
		globals (objects);                                                      /* Sum over processors */

	/* Calculate friction matrix and its inverse of the particle  */

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{	
		sph = &objects[n_sph];
		sph->ztt.xx -= sph->c.x*sph->c.x/sph->sum;
		sph->ztt.yy -= sph->c.y*sph->c.y/sph->sum;
		sph->ztt.zz -= sph->c.z*sph->c.z/sph->sum;
		sph->ztt.yz -= sph->c.y*sph->c.z/sph->sum;
		sph->ztt.zx -= sph->c.z*sph->c.x/sph->sum;
		sph->ztt.xy -= sph->c.x*sph->c.y/sph->sum;
		sph->zrr.xx -= sph->rxc.x*sph->rxc.x/sph->sum;
		sph->zrr.yy -= sph->rxc.y*sph->rxc.y/sph->sum;
		sph->zrr.zz -= sph->rxc.z*sph->rxc.z/sph->sum;
		sph->zrr.yz -= sph->rxc.y*sph->rxc.z/sph->sum;
		sph->zrr.zx -= sph->rxc.z*sph->rxc.x/sph->sum;
		sph->zrr.xy -= sph->rxc.x*sph->rxc.y/sph->sum;
		sph->ztr.xx -= sph->c.x*sph->rxc.x/sph->sum;
		sph->ztr.yy -= sph->c.y*sph->rxc.y/sph->sum;
		sph->ztr.zz -= sph->c.z*sph->rxc.z/sph->sum;
		sph->ztr.yz -= sph->c.y*sph->rxc.z/sph->sum;
		sph->ztr.zx -= sph->c.z*sph->rxc.x/sph->sum;
		sph->ztr.xy -= sph->c.x*sph->rxc.y/sph->sum;
		sph->ztr.zy -= sph->c.z*sph->rxc.y/sph->sum;
		sph->ztr.xz -= sph->c.x*sph->rxc.z/sph->sum;
		sph->ztr.yx -= sph->c.y*sph->rxc.x/sph->sum;

		if (sph->mass_flag)  matrix_inv (sph);                              /* m-matrices from inverse of z's */
	}
}


/*  MATRIX_INV: Inverse of Boundary node moments */

void matrix_inv (struct object *sph)
{
	double z[6][12];
	int    i, j;

	z[0][0] = sph->ztt.xx + sph->mass;
	z[0][1] = sph->ztt.xy;
	z[0][2] = sph->ztt.zx;
	z[1][0] = sph->ztt.xy;
	z[1][1] = sph->ztt.yy + sph->mass;
	z[1][2] = sph->ztt.yz;
	z[2][0] = sph->ztt.zx;
	z[2][1] = sph->ztt.yz;
	z[2][2] = sph->ztt.zz + sph->mass;

	z[3][3] = sph->zrr.xx + sph->inertia;
	z[3][4] = sph->zrr.xy;
	z[3][5] = sph->zrr.zx;
	z[4][3] = sph->zrr.xy;
	z[4][4] = sph->zrr.yy + sph->inertia;
	z[4][5] = sph->zrr.yz;
	z[5][3] = sph->zrr.zx;
	z[5][4] = sph->zrr.yz;
	z[5][5] = sph->zrr.zz + sph->inertia;

	z[0][3] = sph->ztr.xx;
	z[0][4] = sph->ztr.xy;
	z[0][5] = sph->ztr.xz;
	z[1][3] = sph->ztr.yx;
	z[1][4] = sph->ztr.yy;
	z[1][5] = sph->ztr.yz;
	z[2][3] = sph->ztr.zx;
	z[2][4] = sph->ztr.zy;
	z[2][5] = sph->ztr.zz;

	z[3][0] = sph->ztr.xx;
	z[3][1] = sph->ztr.yx;
	z[3][2] = sph->ztr.zx;
	z[4][0] = sph->ztr.xy;
	z[4][1] = sph->ztr.yy;
	z[4][2] = sph->ztr.zy;
	z[5][0] = sph->ztr.xz;
	z[5][1] = sph->ztr.yz;
	z[5][2] = sph->ztr.zz;

	for (i = 0; i <  6; i++)
		for (j = 6; j < 12; j++)
			z[i][j] = 0.0;
	for (i = 0; i <  6; i++)
		z[i][6+i] = 1.0;

	gauss_jordan (z);

	sph->mtt.xx = z[0][ 6];
	sph->mtt.yy = z[1][ 7];
	sph->mtt.zz = z[2][ 8];
	sph->mtt.yz = z[1][ 8];
	sph->mtt.zx = z[2][ 6];
	sph->mtt.xy = z[0][ 7];

	sph->mrr.xx = z[3][ 9];
	sph->mrr.yy = z[4][10];
	sph->mrr.zz = z[5][11];
	sph->mrr.yz = z[4][11];
	sph->mrr.zx = z[5][ 9];
	sph->mrr.xy = z[3][10];

	sph->mtr.xx = z[0][ 9];
	sph->mtr.xy = z[0][10];
	sph->mtr.xz = z[0][11];
	sph->mtr.yx = z[1][ 9];
	sph->mtr.yy = z[1][10];
	sph->mtr.yz = z[1][11];
	sph->mtr.zx = z[2][ 9];
	sph->mtr.zy = z[2][10];
	sph->mtr.zz = z[2][11];
}


/*  GAUSS-JORDAN: Matrix inversion by Gauss-Jordan  */

void gauss_jordan (double z[6][12])
{
	double ratio;
	int    i, j, k;

	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
		{
			if (i == j)  continue;
			ratio = z[j][i]/z[i][i];
			for (k = 0; k < 12; k++)
				z[j][k] -= ratio*z[i][k];
		}
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
			z[i][j+6] /= z[i][i];
}


/*  BNODES_ADD: Create fluid in old bounday node regions  */

void bnodes_add (struct object *objects, Float ***velcs_df)
{
	struct object *sph;
	double rx, ry, rz, ux, uy, uz, wx, wy, wz, dx, dy, dz;
	double rbx, rby, rbz, ubx, uby, ubz, fbx, fby, fbz;
	int    x1, y1, z1, l1, xy1;
	int    node, flag, n_sph, n_bnode;

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{	
		if (objects[n_sph].num_bnode == 0)  continue;	
		sph = &objects[n_sph];
		rx  = sph->r_map.x;
		ry  = sph->r_map.y;
		rz  = sph->r_map.z;
		dx  = sph->r.x - rx;                                                /* New and old coord differences */
		dy  = sph->r.y - ry;
		dz  = sph->r.z - rz;
		ux  = sph->u.x;
		uy  = sph->u.y;
		uz  = sph->u.z;
		wx  = sph->w.x;
		wy  = sph->w.y;
		wz  = sph->w.z;

		for (n_bnode = 0; n_bnode < sph->num_bnode; n_bnode++)
		{
			node = sph->nodes[n_bnode];                                     /* Boundary nodes */
			x1   = node/((max_y+2)*(max_z+2)) - 1;                                  /* get node index for the coordinates */
			y1   = node%((max_y+2)*(max_z+2))/(max_z+2);
			z1   = node%(max_z+2);
			xy1  = (x1+1)*(max_y+2) + y1;
			if (range (x1, 0, num_x))  flag = 1;                            /* Node in processor domain */
			else                       flag = 0;
			rbx  = x1 + 0.5 - rx + x_min;                                   /* Coords (incl. box start) */
			rby  = y1 + 0.5 - ry;
			rbz  = z1 + 0.5 - rz;
			rbx  = n_image(rbx, max_x);                                     /* Periodic BC's */
			rby  = n_image(rby, max_y);
			rbz  = n_image(rbz, max_z);

			/* Particle boundary condition of mass transfer across the solid-fluid boundary at the outside node */
			if (sph->func_ptr(sph, rbx-dx, rby-dy, rbz-dz))                 /* Only nodes outside new surface */
			{
				/* u_boundary = u + w cross rb */
				ubx  = ux + wy*rbz - wz*rby;
				uby  = uy + wz*rbx - wx*rbz;
				ubz  = uz + wx*rby - wy*rbx;
				fbx  = Rho_Fl*ubx*flag;
				fby  = Rho_Fl*uby*flag;
				fbz  = Rho_Fl*ubz*flag;

				/* df = a * (ub dot c) */
				for (l1 = 0; l1 < Num_Dir; l1++)
					velcs_df[xy1][l1][z1] = fac[l1]*(ubx*c_x[l1] + uby*c_y[l1] + ubz*c_z[l1])/CS2;

				sph->f.x -= fbx;
				sph->f.y -= fby;
				sph->f.z -= fbz;
				sph->t.x -= rby*fbz - rbz*fby;
				sph->t.y -= rbz*fbx - rbx*fbz;
				sph->t.z -= rbx*fby - rby*fbx;
			}
		}
	}
}


/*  BNODES_DEL: Delete fluid in new bounday node regions  */

void bnodes_del (struct object *objects, Float ***velcs_df)
{
	struct object *sph;
	double rx, ry, rz;
	double rbx, rby, rbz, fbx, fby, fbz, fb, mass;
	int    x1, y1, z1, l1, xy1;
	int    node, flag, n_sph, n_bnode;

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{	
		if (objects[n_sph].num_bnode == 0)  continue;	
		sph  = &objects[n_sph];
		rx   = sph->r_map.x;
		ry   = sph->r_map.y;
		rz   = sph->r_map.z;
		mass = 0.0;

		for (n_bnode = 0; n_bnode < sph->num_bnode; n_bnode++)
		{
			node = sph->nodes[n_bnode];                                     /* Boundary nodes */
			x1   = node/((max_y+2)*(max_z+2)) - 1;                                  /* get node index from coordinates */
			y1   = node%((max_y+2)*(max_z+2))/(max_z+2);
			z1   = node%(max_z+2);
			xy1  = (x1+1)*(max_y+2) + y1;
			if (range (x1, 0, num_x))  flag = 1;                            /* Node in processor domain */
			else                       flag = 0;
			rbx  = x1 + 0.5 - rx + x_min;                                   /* Coords (incl. box start) */
			rby  = y1 + 0.5 - ry;
			rbz  = z1 + 0.5 - rz;
			rbx  = n_image(rbx, max_x);                                     /* Periodic BC's */
			rby  = n_image(rby, max_y);
			rbz  = n_image(rbz, max_z);
			fbx  = 0.0;
			fby  = 0.0;
			fbz  = 0.0;

			for (l1 = 0; l1 < Num_Dir; l1++)
			{
				fb    = velcs_df[xy1][l1][z1]*flag;
				mass += fb;
				fbx  += fb*c_x[l1];
				fby  += fb*c_y[l1];
				fbz  += fb*c_z[l1];
				velcs_df[xy1][l1][z1] = 0.0;
			}

			sph->f.x += fbx;
			sph->f.y += fby;
			sph->f.z += fbz;
			sph->t.x += rby*fbz - rbz*fby;
			sph->t.y += rbz*fbx - rbx*fbz;
			sph->t.z += rbx*fby - rby*fbx;
		}
		sph->fluid += mass;
	}
}
