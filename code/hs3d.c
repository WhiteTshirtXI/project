/*  HS3D: Hard sphere MD code: updates system for dt  */
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

extern int  num_proc, n_proc;
extern int  max_x, max_y, max_z;
extern int  num_sph, num_obj;

void hs3d (struct object *objects, double dt)
{
  struct list *list;
  double  x12, y12, z12, r12, sum_ke1, sum_ke2;
  double  *t_last, t_next, t_now;
  double  m1, m2, m, t1, t2, w1, w2, sig;
  double  ux1, uy1, uz1, ux2, uy2, uz2;
  double  ur, urx, ury, urz, prx, pry, prz;
  int     n1_sph, n2_sph, f1, f2;
  int     num_sph_coll, num_wall_coll;
  int     first, last, now, next, n_cycle;
	
	   
	
  dt /= (double) Num_Hs3d_Step;
  for (n_cycle = 0; n_cycle < Num_Hs3d_Step; n_cycle++)
    {
      t_last = (double *) calloc (num_sph, sizeof(*t_last));
      list   = (struct list *) calloc (MAX_C, sizeof(struct list));
      if (t_last == 0)  fatal_err ("cannot allocate t_last", num_sph);
      if (list == 0)    fatal_err ("cannot allocate list", MAX_C);
	   
      /*  Create initial collision list    */
	   
      last = 0;   t_now = 0.0;   sum_ke1 = 0.0;
      for (n1_sph = 0; n1_sph < num_sph; n1_sph++)
	{
	  m1  = objects[n1_sph].mass*objects[n1_sph].mass_flag;
	  ux1 = objects[n1_sph].u.x;
	  uy1 = objects[n1_sph].u.y;
	  uz1 = objects[n1_sph].u.z;
	  sum_ke1 += 0.5*m1*(ux1*ux1 + uy1*uy1 + uz1*uz1);
	  last = coll (objects, list, t_last, t_now, dt, n1_sph, 0, n1_sph, last);
	}

      /*  Do all collisions for next dt  */
	   
      first = 0;   t_next = 0.0;   num_sph_coll = 0;  num_wall_coll = 0;
      sum_ke2 = 0.0;
      while (t_next < dt)
	{
	  t_next = dt + 0.1;                                    /* Find next collision */
	  now    = first;
	  while (now < last)
	    {
	      if (list[now].t_coll < t_next)
		{
		  n1_sph = list[now].coll_1;
		  n2_sph = list[now].coll_2;
		  t_next = list[now].t_coll;
		}
	      now = list[now].next;
	    }
		
	  if (t_next < dt)
	    {
	      t_now = t_next;
		
	      m1  = objects[n1_sph].mass;
	      f1  =(objects[n1_sph].move_flag+1)/2;
	      t1  = t_now - t_last[n1_sph];
	      ux1 = objects[n1_sph].u.x*f1;
	      uy1 = objects[n1_sph].u.y*f1;
	      uz1 = objects[n1_sph].u.z*f1;
	      objects[n1_sph].r.x += ux1*t1;
	      objects[n1_sph].r.y += uy1*t1;
	      objects[n1_sph].r.z += uz1*t1;
	      objects[n1_sph].p_str.xx -= m1*ux1*ux1*t1;
	      objects[n1_sph].p_str.yy -= m1*uy1*uy1*t1;
	      objects[n1_sph].p_str.zz -= m1*uz1*uz1*t1;
	      objects[n1_sph].p_str.yz -= m1*uy1*uz1*t1;
	      objects[n1_sph].p_str.zx -= m1*uz1*ux1*t1;
	      objects[n1_sph].p_str.xy -= m1*ux1*uy1*t1;
	      t_last[n1_sph] = t_now;
			
	      if (n2_sph < num_sph)	                          /* Update particle-particle velocities */
		{
		  m2  = objects[n2_sph].mass;
		  f2  =(objects[n2_sph].move_flag+1)/2;
		  t2  = t_now - t_last[n2_sph];
		  ux2 = objects[n2_sph].u.x*f2;
		  uy2 = objects[n2_sph].u.y*f2;
		  uz2 = objects[n2_sph].u.z*f2;
		  objects[n2_sph].r.x += ux2*t2;
		  objects[n2_sph].r.y += uy2*t2;
		  objects[n2_sph].r.z += uz2*t2;
		  objects[n2_sph].p_str.xx -= m2*ux2*ux2*t2;
		  objects[n2_sph].p_str.yy -= m2*uy2*uy2*t2;
		  objects[n2_sph].p_str.zz -= m2*uz2*uz2*t2;
		  objects[n2_sph].p_str.yz -= m2*uy2*uz2*t2;
		  objects[n2_sph].p_str.zx -= m2*uz2*ux2*t2;
		  objects[n2_sph].p_str.xy -= m2*ux2*uy2*t2;
		  t_last[n2_sph] = t_now;
					
		  x12 = n_image(objects[n1_sph].r.x - objects[n2_sph].r.x, max_x);
		  y12 = n_image(objects[n1_sph].r.y - objects[n2_sph].r.y, max_y);
		  z12 = n_image(objects[n1_sph].r.z - objects[n2_sph].r.z, max_z);
		  r12 = x12*x12 + y12*y12 + z12*z12;
		  sig = objects[n1_sph].r_a + objects[n2_sph].r_a;
		  if (fabs(r12-sig*sig) > Tol)  warning ("illegal collision");
		  urx = ux1 - ux2;
		  ury = uy1 - uy2;
		  urz = uz1 - uz2;
		  ur  = (urx*x12 + ury*y12 + urz*z12)/sig;

		  if (objects[n1_sph].move_flag == 2)  f1 = 0;  /* Wall particles unaffected by collisions */
		  if (objects[n2_sph].move_flag == 2)  f2 = 0;

		  m   = 1.0/(f1/m1 + f2/m2);                    /* Reduced mass */
		  prx = m*ur*(1.0 + Elas_C)*x12/sig;            /* Inelastic collision */
		  pry = m*ur*(1.0 + Elas_C)*y12/sig;
		  prz = m*ur*(1.0 + Elas_C)*z12/sig;
		  objects[n1_sph].u.x -= prx*f1/m1;
		  objects[n1_sph].u.y -= pry*f1/m1;
		  objects[n1_sph].u.z -= prz*f1/m1;
		  objects[n2_sph].u.x += prx*f2/m2;
		  objects[n2_sph].u.y += pry*f2/m2;
		  objects[n2_sph].u.z += prz*f2/m2;
		  objects[n1_sph].p.x -= prx;
		  objects[n1_sph].p.y -= pry;
		  objects[n1_sph].p.z -= prz;
		  objects[n2_sph].p.x += prx;
		  objects[n2_sph].p.y += pry;
		  objects[n2_sph].p.z += prz;

		  w1 = (f1/m1)/(f1/m1+f2/m2);
		  w2 = (f2/m2)/(f1/m1+f2/m2);
		  objects[n1_sph].pc_str.xx += prx*x12*w1;
		  objects[n1_sph].pc_str.yy += pry*y12*w1;
		  objects[n1_sph].pc_str.zz += prz*z12*w1;
		  objects[n1_sph].pc_str.yz += pry*z12*w1;
		  objects[n1_sph].pc_str.zx += prz*x12*w1;
		  objects[n1_sph].pc_str.xy += prx*y12*w1;
		  objects[n2_sph].pc_str.xx += prx*x12*w2;
		  objects[n2_sph].pc_str.yy += pry*y12*w2;
		  objects[n2_sph].pc_str.zz += prz*z12*w2;
		  objects[n2_sph].pc_str.yz += pry*z12*w2;
		  objects[n2_sph].pc_str.zx += prz*x12*w2;
		  objects[n2_sph].pc_str.xy += prx*y12*w2;
		  objects[n1_sph].e_diss.c  += 0.5*m*(1.0 - Elas_C*Elas_C)*ur*ur*w1;
		  objects[n2_sph].e_diss.c  += 0.5*m*(1.0 - Elas_C*Elas_C)*ur*ur*w2;
		  if (f1 == 0)  sum_ke2 -= ux1*prx + uy1*pry + uz1*prz;
		  if (f2 == 0)  sum_ke2 += ux2*prx + uy2*pry + uz2*prz;
		  objects[n1_sph].num_coll++;
		  objects[n2_sph].num_coll++;
		  num_sph_coll++;
		}
	      else
		{
		  x12 =(objects[n1_sph].r.x - objects[n2_sph].r.x)*objects[n2_sph].e.x*objects[n2_sph].e.x;
		  y12 =(objects[n1_sph].r.y - objects[n2_sph].r.y)*objects[n2_sph].e.y*objects[n2_sph].e.y;
		  z12 =(objects[n1_sph].r.z - objects[n2_sph].r.z)*objects[n2_sph].e.z*objects[n2_sph].e.z;
		  sig = objects[n1_sph].r_a;
		  if (x12+y12+z12-sig > Tol)  warning ("illegal wall collision");
		  urx = ux1;
		  ury = uy1;
		  urz = uz1;
		  ur  = (urx*x12 + ury*y12 + urz*z12)/sig;
		  m   = m1;                                     /* Reduced mass */
		  prx = m*ur*(1.0 + Elas_C)*x12/sig;            /* Inelastic collision */
		  pry = m*ur*(1.0 + Elas_C)*y12/sig;
		  prz = m*ur*(1.0 + Elas_C)*z12/sig;
		  objects[n1_sph].u.x -= prx/m1;
		  objects[n1_sph].u.y -= pry/m1;
		  objects[n1_sph].u.z -= prz/m1;
		  objects[n1_sph].p.x -= prx;
		  objects[n1_sph].p.y -= pry;
		  objects[n1_sph].p.z -= prz;
		  objects[n2_sph].p.x += prx;
		  objects[n2_sph].p.y += pry;
		  objects[n2_sph].p.z += prz;
		  objects[n1_sph].pc_str.xx += prx*x12;
		  objects[n1_sph].pc_str.yy += pry*y12;
		  objects[n1_sph].pc_str.zz += prz*z12;
		  objects[n1_sph].e_diss.c  += 0.5*m*(1.0 - Elas_C*Elas_C)*ur*ur;
		  objects[n1_sph].num_coll++;
					
		  num_wall_coll++;
		}
			
	      /*	update collision time table  */
			
	      if (n2_sph < num_sph)
		{
		  now  = first;
		  next = list[now].next;
		  while (next < last)
		    {
		      if ((list[next].coll_1 == n1_sph) || (list[next].coll_2 == n1_sph)
			  ||	(list[next].coll_1 == n2_sph) || (list[next].coll_2 == n2_sph))  list[now].next  = list[next].next;
		      else  now  = next;
		      next = list[next].next;
		    }
		  if ((list[first].coll_1 == n1_sph) || (list[first].coll_2 == n1_sph)
		      ||	(list[first].coll_1 == n2_sph) || (list[first].coll_2 == n2_sph))  first = list[first].next;
		  last = coll (objects, list, t_last, t_now, dt, n1_sph, 0, num_sph, last);
		  last = coll (objects, list, t_last, t_now, dt, n2_sph, 0, num_sph, last);
		}
	      else 
		{
		  now  = first;
		  next = list[now].next;
		  while (next < last)
		    {
		      next = list[now].next;
		      if ((list[next].coll_1 == n1_sph) || (list[next].coll_2 == n1_sph))  list[now].next  = list[next].next;
		      else  now  = next;
		      next = list[next].next;
		    }
		  if ((list[first].coll_1 == n1_sph) || (list[first].coll_2 == n1_sph))  first = list[first].next;
		  last = coll (objects, list, t_last, t_now, dt, n1_sph, 0, num_sph, last);
		}
	    }
	}

      /*	 return new coordinates  */

      for (n1_sph = 0; n1_sph < num_sph; n1_sph++)
	{
	  m1  = objects[n1_sph].mass*objects[n1_sph].mass_flag;
	  ux1 = objects[n1_sph].u.x;
	  uy1 = objects[n1_sph].u.y;
	  uz1 = objects[n1_sph].u.z;
	  sum_ke2 += 0.5*m1*(ux1*ux1 + uy1*uy1 + uz1*uz1) + objects[n1_sph].e_diss.c;
	  if (ux1*ux1+uy1*uy1+uz1*uz1 > 1.0) fatal_err ("Particle moving too fast", n1_sph);

	  m1  = objects[n1_sph].mass;
	  f1  =(objects[n1_sph].move_flag+1)/2;
	  t1  = dt-t_last[n1_sph];
	  ux1 *= f1;
	  uy1 *= f1;
	  uz1 *= f1;
	  objects[n1_sph].r.x += ux1*t1;
	  objects[n1_sph].r.y += uy1*t1;
	  objects[n1_sph].r.z += uz1*t1;
	  objects[n1_sph].p_str.xx -= m1*ux1*ux1*t1;
	  objects[n1_sph].p_str.yy -= m1*uy1*uy1*t1;
	  objects[n1_sph].p_str.zz -= m1*uz1*uz1*t1;
	  objects[n1_sph].p_str.yz -= m1*uy1*uz1*t1;
	  objects[n1_sph].p_str.zx -= m1*uz1*ux1*t1;
	  objects[n1_sph].p_str.xy -= m1*ux1*uy1*t1;
	}

      free (t_last);
      free (list);

      if (fabs(sum_ke2-sum_ke1)/(sum_ke1 + Tol) > Tol)
	{
	  fprintf (stdout, "Initial KE = % .5e: Final KE = % .5e\n", sum_ke1, sum_ke2);
	  fatal_err ("Energy conservation error", -1);
	}
      /*		if (n_proc == 0)  fprintf (stdout, " %d particle-particle collision(s): %d particle-wall collision(s)\n", num_sph_coll, num_wall_coll);  */
    }

	

}


int coll (struct object *objects, struct list *list, double *t_last, double t_now, double dt, int n1_sph, int n1, int n2, int last)
{
	double  x1, y1, z1, x2, y2, z2, x12, y12, z12, r12, sig;
	double  ux1, uy1, uz1, ux2, uy2, uz2;
	double  urx, ury, urz, ur, u12, ur2, t12;
	int     n2_sph, n2_obj, f1, f2, n;


	f1  =(objects[n1_sph].move_flag+1)/2;
	ux1 = objects[n1_sph].u.x*f1;
	uy1 = objects[n1_sph].u.y*f1;
	uz1 = objects[n1_sph].u.z*f1;
	x1  = objects[n1_sph].r.x + ux1*(t_now-t_last[n1_sph]);
	y1  = objects[n1_sph].r.y + uy1*(t_now-t_last[n1_sph]);
	z1  = objects[n1_sph].r.z + uz1*(t_now-t_last[n1_sph]);

	for (n = 1; n <= objects[n1_sph].list[0]; n++)
	{	
		n2_sph = objects[n1_sph].list[n];
		if (range(n2_sph, n1, n2))
		{
			f2  =(objects[n2_sph].move_flag+1)/2;
			ux2 = objects[n2_sph].u.x*f2;
			uy2 = objects[n2_sph].u.y*f2;
			uz2 = objects[n2_sph].u.z*f2;
			x2  = objects[n2_sph].r.x + ux2*(t_now-t_last[n2_sph]);
			y2  = objects[n2_sph].r.y + uy2*(t_now-t_last[n2_sph]);
			z2  = objects[n2_sph].r.z + uz2*(t_now-t_last[n2_sph]);
			x12 = n_image(x1 - x2, max_x);
			y12 = n_image(y1 - y2, max_y);
			z12 = n_image(z1 - z2, max_z);
			urx = ux1 - ux2;
			ury = uy1 - uy2;
			urz = uz1 - uz2;
			sig = objects[n1_sph].r_a + objects[n2_sph].r_a;
			r12 = x12*x12 + y12*y12 + z12*z12 - sig*sig;
			if (r12 < -Tol)  warning ("particle-particle overlap");
			
			ur  = urx*x12 + ury*y12 + urz*z12;
			if (ur < 0.0)
			{
				u12 = urx*urx + ury*ury + urz*urz;
				ur2 = u12*r12;
				if (ur*ur >= ur2)
				{
					t12 = - ur - sqrt(ur*ur - ur2);
					t12 = t12/u12 + t_now;
					if (t12 < dt)
					{
						list[last].t_coll = t12;
						list[last].coll_1 = n1_sph;
						list[last].coll_2 = n2_sph;
						list[last].next   = last + 1;
						last++;
						if (last > MAX_C)  fatal_err("too many collisions: MAX_C exceeded", -1);
					}
				}
			}
		}
	}
	
	
	if (objects[n1_sph].move_flag > 0)
	{
		for (n2_obj = num_sph; n2_obj < num_obj; n2_obj++)
		{
			x12 = (x1 - objects[n2_obj].r.x)*objects[n2_obj].e.x;
			y12 = (y1 - objects[n2_obj].r.y)*objects[n2_obj].e.y;
			z12 = (z1 - objects[n2_obj].r.z)*objects[n2_obj].e.z;

			sig =  objects[n1_sph].r_a;
			r12 = x12 + y12 + z12 - sig;
			if (r12 < -Tol)  warning ("particle-wall overlap");
		
			urx = objects[n1_sph].u.x*objects[n2_obj].e.x;
			ury = objects[n1_sph].u.y*objects[n2_obj].e.y;
			urz = objects[n1_sph].u.z*objects[n2_obj].e.z;
			ur	= urx + ury + urz;
			if (ur < 0.0)
			{
				t12 = -r12/ur + t_now;
				if (t12 < dt)
				{
					list[last].t_coll = t12;
					list[last].coll_1 = n1_sph;
					list[last].coll_2 = n2_obj;
					list[last].next   = last + 1;
					last++;
					if (last > MAX_C)  fatal_err("too many collisions: MAX_C exceeded", -1);
				}
			}
		}
	}
	return (last);
}
