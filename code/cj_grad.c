/*  Conjugate-gradient solver A.x = b  */
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

void cj_grad (double **A, double *x, int size)
{

	double *g, *h, *r, bb, hg, rr, rr_old, alpha, beta;
	int    *ind_i, *ind_j, num_ind, index, n, i, j;

	g  = (double *) calloc (size, sizeof(double));
	h  = (double *) calloc (size, sizeof(double));
	r  = (double *) calloc (size, sizeof(double));
	ind_i = (int *) calloc (size*size, sizeof(int));
	ind_j = (int *) calloc (size*size, sizeof(int));
	if (g == 0) fatal_err("cannot allocate g", size);
	if (h == 0) fatal_err("cannot allocate h", size);
	if (r == 0) fatal_err("cannot allocate r", size);
	if (ind_i == 0) fatal_err("cannot allocate ind_i", size*size);
	if (ind_j == 0) fatal_err("cannot allocate ind_j", size*size);
	index = 0;

	for (i = 0; i < size; i++)
	for (j = 0; j < size; j++)
	{
		if (A[i][j] != 0.0)
		{
			ind_i[index] = i;
			ind_j[index] = j;
			index++;
		}
	}
	num_ind = index;

	for (i = 0, bb = 0.0,  rr = 0.0; i < size; i++)
		{
		g[i] = r[i] = x[i];                               /* g0 = r0 = b */
		bb  += x[i]*x[i];                                 /* b norm */
		rr  += r[i]*r[i];                                 /* r norm */
		x[i] = 0.0;                                       /* Initial x = 0 */
		}

	for (n = 1; n <= size && sqrt(rr/bb) > 0.1*Tol; n++)
	{
		for (i = 0; i < size; i++)
			h[i] = 0.0;
		for (index = 0; index < num_ind; index++)         /* h = A.g */
		{
			i = ind_i[index];
			j = ind_j[index];
			h[i] += A[i][j]*g[j];
		}

		for (i = 0, hg = 0.0; i < size; i++)
			hg += h[i]*g[i];

		alpha = rr/hg;                                    /* r.r/h.g */
		for (i = 0; i < size;  i++)
		{
			x[i] += alpha*g[i];                           /* New x */
			r[i] -= alpha*h[i];                           /* New r */
		}

		rr_old = rr;
		for (i = 0, rr = 0.0; i < size; i++)
			rr += r[i]*r[i];

		beta = rr/rr_old;                                 /* r.r/ro.ro */
		for (i = 0; i < size;  i++)
			g[i] = r[i] + beta*g[i];                      /* New g */
	}
	free (ind_i);
	free (ind_j);
	free (g);
	free (h);
	free (r);
}
