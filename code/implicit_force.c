/*  Implicit_Force:  Calculate particle-fluid force using implicit velocity update */
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

void implicit_force (struct object *objects, struct vector f_ext)
{
  extern int  num_obj;
  struct object  *obj;
  struct smatrix *ztt, *zrr, *mtt, *mrr;
  struct matrix  *ztr, *mtr;
  double ux, uy, uz, wx, wy, wz;
  double fx, fy, fz, tx, ty, tz;
  int    n_obj;

  for (n_obj = 0; n_obj < num_obj; n_obj++)
    {
      obj = &objects[n_obj];
      ux  = obj->u.x;
      uy  = obj->u.y;
      uz  = obj->u.z;
      wx  = obj->w.x;
      wy  = obj->w.y;
      wz  = obj->w.z;
      fx  = obj->f.x;
      fy  = obj->f.y;
      fz  = obj->f.z;
      tx  = obj->t.x;
      ty  = obj->t.y;
      tz  = obj->t.z;
      ztt = &obj->ztt;
      ztr = &obj->ztr;
      zrr = &obj->zrr;
      mtt = &obj->mtt;
      mtr = &obj->mtr;
      mrr = &obj->mrr;

      fx -= ztt->xx*ux + ztt->xy*uy + ztt->zx*uz + ztr->xx*wx + ztr->xy*wy + ztr->xz*wz;    /* Drag force */
      fy -= ztt->xy*ux + ztt->yy*uy + ztt->yz*uz + ztr->yx*wx + ztr->yy*wy + ztr->yz*wz;
      fz -= ztt->zx*ux + ztt->yz*uy + ztt->zz*uz + ztr->zx*wx + ztr->zy*wy + ztr->zz*wz;
      tx -= ztr->xx*ux + ztr->yx*uy + ztr->zx*uz + zrr->xx*wx + zrr->xy*wy + zrr->zx*wz;
      ty -= ztr->xy*ux + ztr->yy*uy + ztr->zy*uz + zrr->xy*wx + zrr->yy*wy + zrr->yz*wz;
      tz -= ztr->xz*ux + ztr->yz*uy + ztr->zz*uz + zrr->zx*wx + zrr->yz*wy + zrr->zz*wz;

      
      fx += obj->vol*f_ext.x;                                                               /* Buoyancy force */
      fy += obj->vol*f_ext.y;
      fz += obj->vol*f_ext.z;
      
      obj->e_diss.t -= ux*fx + uy*fy + uz*fz;                                               /* Energy dissipation */
      obj->e_diss.r -= wx*tx + wy*ty + wz*tz;

      fx += obj->f_ext.x;                                                                   /* External force */
      fy += obj->f_ext.y;
      fz += obj->f_ext.z;

      /* Implicit update of the velocity */
      /* U(t+1) = U(t) + inv. friction matrix dot F(t+1) */
      /* See Nguyen and Ladd (2002) */

      if (objects[n_obj].mass_flag)
	{
	  ux = mtt->xx*fx + mtt->xy*fy + mtt->zx*fz + mtr->xx*tx + mtr->xy*ty + mtr->xz*tz; /* Velocity change */
	  uy = mtt->xy*fx + mtt->yy*fy + mtt->yz*fz + mtr->yx*tx + mtr->yy*ty + mtr->yz*tz;
	  uz = mtt->zx*fx + mtt->yz*fy + mtt->zz*fz + mtr->zx*tx + mtr->zy*ty + mtr->zz*tz;
	  wx = mtr->xx*fx + mtr->yx*fy + mtr->zx*fz + mrr->xx*tx + mrr->xy*ty + mrr->zx*tz;
	  wy = mtr->xy*fx + mtr->yy*fy + mtr->zy*fz + mrr->xy*tx + mrr->yy*ty + mrr->yz*tz;
	  wz = mtr->xz*fx + mtr->yz*fy + mtr->zz*fz + mrr->zx*tx + mrr->yz*ty + mrr->zz*tz;
	  obj->u.x += ux;
	  obj->u.y += uy;
	  obj->u.z += uz;
	  obj->w.x += wx;
	  obj->w.y += wy;
	  obj->w.z += wz;
	  fx = ux*obj->mass;                                                                /* Implicit force */
	  fy = uy*obj->mass;
	  fz = uz*obj->mass;
	  tx = wx*obj->inertia;
	  ty = wy*obj->inertia;
	  tz = wz*obj->inertia;
	}

      /* linear and angular momentum update */
      obj->p.x += fx;                                                                       /* Increment average force */
      obj->p.y += fy;
      obj->p.z += fz;
      obj->l.x += tx;
      obj->l.y += ty;
      obj->l.z += tz;
      obj->f.x  = 0.0;
      obj->f.y  = 0.0;
      obj->f.z  = 0.0;
      obj->t.x  = 0.0;
      obj->t.y  = 0.0;
      obj->t.z  = 0.0;
    }
}
