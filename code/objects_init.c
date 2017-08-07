/*  OBJECTS_INIT: Default Initializations (overridden in driver)  */
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

#include  "header.h"

/*  SPHEROID_INIT: Spheroidal particle */

void spheroid_init (struct object *sph, double mass_fac, double vel_fac)
{
  double r_a, r_b;	
	
  r_a = sph->r_a;
  r_b = sph->r_b;
	
  sph->r_map.x   = sph->r.x;
  sph->r_map.y   = sph->r.y;
  sph->r_map.z   = sph->r.z;
  sph->size      = max (r_a, r_b) + 1;
  sph->dr2       = 0.0;
	
  switch (sph->shape_flag)
  {
    case 0:
    sph->func_ptr  = sphere;                              /* pointer to function for object mapping */
    sph->max_bnode = 4*Pi*sqrt(2.0)*(r_a+1.0)*(r_a+1.0);  /* Max # of boundary nodes */
    sph->vol       = 4.0*Pi*r_a*r_a*r_a/3.0;          /* Volume fluctuates slightly; Initial volume */
    sph->mass      = sph->rho*Rho_Fl*sph->vol;        /* Particle mass */
    sph->inertia   = 0.4*sph->mass*r_a*r_a;
    break;
    
    case 1:
    fatal_err ("this object is not yet implemented", -1);
    /* To be Modified: objects_map, velcs_update, hs3d, lub */
    break;
    
    case 2:
    fatal_err ("this object is not yet implemented", -1);
    /* To be Modified: objects_map, velcs_update, hs3d, lub */
    break;
    
    case 3:
    fatal_err ("this object is not yet implemented", -1);
    /* To be Modified: objects_map, velcs_update, hs3d, lub */
    break;
    
    case 4:
    fatal_err ("this object is not yet implemented", -1);
    /* To be Modified: objects_map, velcs_update, hs3d, lub */
    break;
  }

  sph->mass     *= mass_fac;  /* Scale mass and inertia */
  sph->inertia  *= mass_fac;
  sph->u.x      *= vel_fac;   /* Scale velocities */
  sph->u.y      *= vel_fac;
  sph->u.z      *= vel_fac;
  sph->w.x      *= vel_fac;
  sph->w.y      *= vel_fac;
  sph->w.z      *= vel_fac;
}


/*  WALL: Plane wall  */

void wall_init (struct object *sph, int n0)
{
  extern int    max_x, max_y, max_z;
  extern int    num_x;
	
  sph->move_flag  = 0;
  sph->mass_flag  = 0;
	
  switch (n0)
  {
    case 0:                                              /* Y- wall */
    sph->r.x         = 0;
    sph->r.y         = 0;
    sph->r.z         = 0;
    sph->e.x         = 0;
    sph->e.y         = 1;
    sph->e.z         = 0;
    sph->max_bnode   = 5*num_x*max_z;
    break;
    
    case 1:                                              /* Y+ wall */
    sph->r.x         = 0;
    sph->r.y         = max_y;
    sph->r.z         = 0;
    sph->e.x         = 0;
    sph->e.y         =-1;
    sph->e.z         = 0;
    sph->max_bnode   = 5*num_x*max_z;
    break;
    
    case 2:                                              /* Z- wall */
    sph->r.x         = 0;
    sph->r.y         = 0;
    sph->r.z         = 0;
    sph->e.x         = 0;
    sph->e.y         = 0;
    sph->e.z         = 1;
    sph->max_bnode   = 5*num_x*max_y;
    break;
    
    case 3:                                              /* Z+ wall */
    sph->r.x         = 0;
    sph->r.y         = 0;
    sph->r.z         = max_z;
    sph->e.x         = 0;
    sph->e.y         = 0;
    sph->e.z         =-1;
    sph->max_bnode   = 5*num_x*max_y;
    break;
    
    case 4:                                              /* X- wall */
    sph->r.x         = 0;
    sph->r.y         = 0;
    sph->r.z         = 0;
    sph->e.x         = 1;
    sph->e.y         = 0;
    sph->e.z         = 0;
    sph->max_bnode   = 5*max_y*max_z;
    break;
    
    case 5:                                              /* X+ wall */
    sph->r.x         = max_x;
    sph->r.y         = 0;
    sph->r.z         = 0;
    sph->e.x         =-1;
    sph->e.y         = 0;
    sph->e.z         = 0;
    sph->max_bnode   = 5*max_y*max_z;
    break;
  }
}
