/* CLUSTERS:   Cluster analysis; Contains
				CLUSTER_INDEX  (Identify and index clusters)
				CLUSTER_MAKE   (Make cluster lists)  */

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
extern int  n_proc;
extern int  num_sph;

void cluster_index (struct object *objects, int *num_cluster)
{
	int cluster_value, cluster_label, max_branch;
 	int n1, n2, n1_sph, n2_sph, num_cl;
	int done_flag;

	/* Label each cluster from 0 to whatever */
	
	for (n1_sph = 0, cluster_label = 0; n1_sph < num_sph; n1_sph++)
	{
		if (objects[n1_sph].cl_list[0] > 0 && objects[n1_sph].mass_flag)
		{
			objects[n1_sph].n_cl = cluster_label;
			cluster_label++;
		}	
		else
			objects[n1_sph].n_cl = -1;
  	}

	/* Organize clusters into cluster index */
	
	max_branch  = 0;
	done_flag = 1;
	while (done_flag)
	{	
		done_flag = 0;
		for (n1_sph = 0; n1_sph < num_sph; n1_sph++)
		{
			if (objects[n1_sph].mass_flag)
			{
				for (n2 = 0; n2 < objects[n1_sph].cl_list[0]; n2++)
				{
					n2_sph = objects[n1_sph].cl_list[n2+1];
					if (objects[n2_sph].mass_flag)
					{	
						if (objects[n1_sph].n_cl > objects[n2_sph].n_cl)
						{
							objects[n1_sph].n_cl = objects[n2_sph].n_cl;
							done_flag = 1;
						}
					}
				}
			}
		}
		max_branch++;
	}

	/* Relabel clusters from 0 to num_cluster */
		
	cluster_value = -1;
	num_cl = 0;
	for (n1_sph = 0; n1_sph < num_sph; n1_sph++)
	{	
		if (objects[n1_sph].n_cl > cluster_value)
		{	
			cluster_value = objects[n1_sph].n_cl;
			objects[n1_sph].n_cl = num_cl;
			num_cl++;
		}
	}					
	done_flag = 1;
	while (done_flag)
	{
		done_flag = 0;
		for (n1_sph = 0; n1_sph < num_sph; n1_sph++)
		{
			if (objects[n1_sph].mass_flag)
			{
				for (n1 = 0; n1 < objects[n1_sph].cl_list[0]; n1++)
				{
					n2_sph = objects[n1_sph].cl_list[n1+1];
					if (objects[n2_sph].mass_flag)
					{
						if (objects[n1_sph].n_cl > objects[n2_sph].n_cl)
						{
							objects[n1_sph].n_cl = objects[n2_sph].n_cl;
							done_flag = 1;
						}	
					}
				}	
			}
		}
	}
	*num_cluster = num_cl;
}

void cluster_make (struct object *objects, struct cluster *clusters, int num_cl)
{
	extern int  max_cl_size;
	int         *list0, n_sph, n_cl;

	list0 = (int *) calloc (num_cl, sizeof(int));
	if (list0 == 0) fatal_err("cannot allocate list0", num_cl);

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{
		if (objects[n_sph].n_cl > -1)
		{
			n_cl = objects[n_sph].n_cl;
			list0[n_cl]++;
			if (list0[n_cl] > max_cl_size)  max_cl_size++;
		}
	}

	for (n_cl = 0; n_cl < num_cl; n_cl++)
	{
		clusters[n_cl].list = (int *) calloc (list0[n_cl]+1, sizeof(int));
		if (clusters[n_cl].list == 0) fatal_err("cannot allocate clusters[n_cl].list", list0[n_cl]+1);
	}

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{
		if (objects[n_sph].n_cl > -1)
		{
			n_cl = objects[n_sph].n_cl;
			clusters[n_cl].list[0]++;
			clusters[n_cl].list[clusters[n_cl].list[0]] = n_sph;
		}
	}

	free (list0);
}
