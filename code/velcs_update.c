/*  VELCS_UPDATE: Updates velocities  */
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

void velcs_update (struct object *objects, double dt)
{
  extern int    n_proc;
  extern int    num_sph, num_obj;
  struct object *sph, *obj;
  struct cluster *clusters;
  int    n_sph, n_obj;
  int    num_cl, n_cl;

  lub (objects, dt);                                                               /* Lubrication forces */
  
  cluster_index (objects, &num_cl);                                               /* Identify clusters */
  if(num_cl > 0) {
    clusters = (struct cluster *) calloc (num_cl, sizeof(struct cluster));
    if (clusters == 0)  fatal_err ("cannot allocate clusters", num_cl);

    cluster_make (objects, clusters, num_cl);                                       /* Make cluster lists */
    
    cluster_force (objects, clusters, dt, num_cl, -1);                              /* Subtract singular forces */
  }

  for (n_sph = 0; n_sph < num_sph; n_sph++)
    {
      sph = &objects[n_sph];
      if (sph->mass_flag)
	{
	  sph->u.x += sph->f.x/sph->mass;                                         /* Update velocities */
	  sph->u.y += sph->f.y/sph->mass;
	  sph->u.z += sph->f.z/sph->mass;
	  sph->w.x += sph->t.x/sph->inertia;
	  sph->w.y += sph->t.y/sph->inertia;
	  sph->w.z += sph->t.z/sph->inertia;
	}
    }
	
  for (n_cl = 0; n_cl < num_cl; n_cl++)
    cluster_update (objects, clusters, dt, n_cl);                               /* Implicit velocity update */

  if(n_cl > 0)
    cluster_force (objects, clusters, dt, num_cl, 1);                               /* Add back singular forces */

  for (n_obj = 0; n_obj < num_obj; n_obj++)
    {
      obj = &objects[n_obj];
      obj->p.x += obj->f.x;                                                       /* Increment average forces */
      obj->p.y += obj->f.y;
      obj->p.z += obj->f.z;
      obj->l.x += obj->t.x;
      obj->l.y += obj->t.y;
      obj->l.z += obj->t.z;
      obj->e_diss.t -= obj->u.x*obj->f.x + obj->u.y*obj->f.y + obj->u.z*obj->f.z; /* Lubrication dissipation */
      obj->e_diss.r -= obj->w.x*obj->t.x + obj->w.y*obj->t.y + obj->w.z*obj->t.z;
      obj->f.x  = 0.0;                                                            /* Zero forces */
      obj->f.y  = 0.0;
      obj->f.z  = 0.0;
      obj->t.x  = 0.0;
      obj->t.y  = 0.0;
      obj->t.z  = 0.0;
    }

  if(n_cl > 0) {
    for (n_cl = 0; n_cl < num_cl; n_cl++)
      free (clusters[n_cl].list);
    free (clusters);
  }
}
