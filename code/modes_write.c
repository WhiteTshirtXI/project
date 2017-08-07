/*  MODES_WRITE: Writes out mass and momentum densities  */
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

void modes_write (Float **velcz_df, int *node_map, struct vector f_ext, FILE *file_ptr)
{
	extern int    max_z;
	int    z, n_dir;
	double m0, m1, m2, m3;
	
	
	/*  Compute normal modes: m0 = mass; m1-m3 = momenta  */
	
	for (z = 0; z < max_z; z++)
	{
		m0 = m1 = m2 = m3 = 0.0;
		if (node_map[z]==0)               /* Fluid nodes */
		{
			for (n_dir = 0; n_dir < Num_Dir; n_dir++) /* Generic LBE */
			{
				m0 += velcz_df[n_dir][z];
				m1 += velcz_df[n_dir][z]*c_x[n_dir];
				m2 += velcz_df[n_dir][z]*c_y[n_dir];
				m3 += velcz_df[n_dir][z]*c_z[n_dir];
			}
			m0 += Rho_Fl;
			m1 -= 0.5*f_ext.x;                        /* Subtract 1/2 total force */
			m2 -= 0.5*f_ext.y;
			m3 -= 0.5*f_ext.z;
		}
		if (file_ptr != 0)               /* Print mass and momentum modes */
			fprintf(file_ptr, " % .5e % .5e % .5e % .5e\n", m0,m1,m2,m3);
	}
}
