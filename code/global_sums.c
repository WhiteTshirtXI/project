/*  GLOBAL_SUMS:  Global sums: Contains
                  GLOBALS    Sum particle properties
                  FORCE_SUM  Global sum of forces and torques
                  MULTI_SUM  Global sum using multiple root nodes  */
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

extern int    num_sph, num_obj;
extern int    num_proc, n_proc;


/*  GLOBALS: Sum particle volume and boundary-node moments  */

void globals (struct object *objects)
{
  double *msg_buf;
  int     num_buf, num_blk, n_buf, n_sph, buf_size, n;

  buf_size =  num_proc*MAX_B;                                         /* Target buffer size */
  msg_buf  = (double *) calloc (buf_size, sizeof(double));
  num_blk  =  buf_size/30;
  num_buf  = (num_sph-1)/num_blk + 1;

  for (n_buf = 0; n_buf < num_buf; n_buf++)
    {
      buf_size = min(num_sph-n_buf*num_blk, num_blk);
      for (n = 0; n < buf_size; n++)
	{
	  n_sph = n_buf*num_blk + n;
	  msg_buf[n*30   ] = objects[n_sph].vol;
	  msg_buf[n*30+ 1] = objects[n_sph].fluid;
	  msg_buf[n*30+ 2] = objects[n_sph].sum;
	  msg_buf[n*30+ 3] = objects[n_sph].c.x;
	  msg_buf[n*30+ 4] = objects[n_sph].c.y;
	  msg_buf[n*30+ 5] = objects[n_sph].c.z;
	  msg_buf[n*30+ 6] = objects[n_sph].rxc.x;
	  msg_buf[n*30+ 7] = objects[n_sph].rxc.y;
	  msg_buf[n*30+ 8] = objects[n_sph].rxc.z;
	  msg_buf[n*30+ 9] = objects[n_sph].ztt.xx;
	  msg_buf[n*30+10] = objects[n_sph].ztt.yy;
	  msg_buf[n*30+11] = objects[n_sph].ztt.zz;
	  msg_buf[n*30+12] = objects[n_sph].ztt.yz;
	  msg_buf[n*30+13] = objects[n_sph].ztt.zx;
	  msg_buf[n*30+14] = objects[n_sph].ztt.xy;
	  msg_buf[n*30+15] = objects[n_sph].zrr.xx;
	  msg_buf[n*30+16] = objects[n_sph].zrr.yy;
	  msg_buf[n*30+17] = objects[n_sph].zrr.zz;
	  msg_buf[n*30+18] = objects[n_sph].zrr.yz;
	  msg_buf[n*30+19] = objects[n_sph].zrr.zx;
	  msg_buf[n*30+20] = objects[n_sph].zrr.xy;
	  msg_buf[n*30+21] = objects[n_sph].ztr.xx;
	  msg_buf[n*30+22] = objects[n_sph].ztr.yy;
	  msg_buf[n*30+23] = objects[n_sph].ztr.zz;
	  msg_buf[n*30+24] = objects[n_sph].ztr.yz;
	  msg_buf[n*30+25] = objects[n_sph].ztr.zx;
	  msg_buf[n*30+26] = objects[n_sph].ztr.xy;
	  msg_buf[n*30+27] = objects[n_sph].ztr.zy;
	  msg_buf[n*30+28] = objects[n_sph].ztr.xz;
	  msg_buf[n*30+29] = objects[n_sph].ztr.yx;
	}

      multi_sum (msg_buf, buf_size*30);

      for (n = 0; n < buf_size; n++)
	{
	  n_sph = n_buf*num_blk + n;
	  objects[n_sph].vol    = msg_buf[n*30   ];
	  objects[n_sph].fluid  = msg_buf[n*30+ 1];
	  objects[n_sph].sum    = msg_buf[n*30+ 2];
	  objects[n_sph].c.x    = msg_buf[n*30+ 3];
	  objects[n_sph].c.y    = msg_buf[n*30+ 4];
	  objects[n_sph].c.z    = msg_buf[n*30+ 5];
	  objects[n_sph].rxc.x  = msg_buf[n*30+ 6];
	  objects[n_sph].rxc.y  = msg_buf[n*30+ 7];
	  objects[n_sph].rxc.z  = msg_buf[n*30+ 8];
	  objects[n_sph].ztt.xx = msg_buf[n*30+ 9];
	  objects[n_sph].ztt.yy = msg_buf[n*30+10];
	  objects[n_sph].ztt.zz = msg_buf[n*30+11];
	  objects[n_sph].ztt.yz = msg_buf[n*30+12];
	  objects[n_sph].ztt.zx = msg_buf[n*30+13];
	  objects[n_sph].ztt.xy = msg_buf[n*30+14];
	  objects[n_sph].zrr.xx = msg_buf[n*30+15];
	  objects[n_sph].zrr.yy = msg_buf[n*30+16];
	  objects[n_sph].zrr.zz = msg_buf[n*30+17];
	  objects[n_sph].zrr.yz = msg_buf[n*30+18];
	  objects[n_sph].zrr.zx = msg_buf[n*30+19];
	  objects[n_sph].zrr.xy = msg_buf[n*30+20];
	  objects[n_sph].ztr.xx = msg_buf[n*30+21];
	  objects[n_sph].ztr.yy = msg_buf[n*30+22];
	  objects[n_sph].ztr.zz = msg_buf[n*30+23];
	  objects[n_sph].ztr.yz = msg_buf[n*30+24];
	  objects[n_sph].ztr.zx = msg_buf[n*30+25];
	  objects[n_sph].ztr.xy = msg_buf[n*30+26];
	  objects[n_sph].ztr.zy = msg_buf[n*30+27];
	  objects[n_sph].ztr.xz = msg_buf[n*30+28];
	  objects[n_sph].ztr.yx = msg_buf[n*30+29];
	}
    }
  free (msg_buf);
}


/*  FORCE_SUM: Sums force contributions over all processors  */

void force_sum (struct object *objects)
{
  double *msg_buf;
  int     num_buf, num_blk, n_buf, n_obj, buf_size, n;

  buf_size =  num_proc*MAX_B;                                         /* Target buffer size */
  msg_buf  = (double *) calloc (buf_size, sizeof(double));
  num_blk  =  buf_size/6;
  num_buf = (num_obj-1)/num_blk + 1;

  for (n_buf = 0; n_buf < num_buf; n_buf++)
    {
      buf_size = min(num_obj-n_buf*num_blk, num_blk);
      for (n = 0; n < buf_size; n++)
	{
	  n_obj = n_buf*num_blk + n;
	  msg_buf[n*6  ] = objects[n_obj].f.x;
	  msg_buf[n*6+1] = objects[n_obj].f.y;
	  msg_buf[n*6+2] = objects[n_obj].f.z;
	  msg_buf[n*6+3] = objects[n_obj].t.x;
	  msg_buf[n*6+4] = objects[n_obj].t.y;
	  msg_buf[n*6+5] = objects[n_obj].t.z;
	}
		
      multi_sum (msg_buf, buf_size*6);
		
      for (n = 0; n < buf_size; n++)
	{
	  n_obj = n_buf*num_blk + n;
	  objects[n_obj].f.x = msg_buf[n*6  ];
	  objects[n_obj].f.y = msg_buf[n*6+1];
	  objects[n_obj].f.z = msg_buf[n*6+2];
	  objects[n_obj].t.x = msg_buf[n*6+3];
	  objects[n_obj].t.y = msg_buf[n*6+4];
	  objects[n_obj].t.z = msg_buf[n*6+5];
	}
    }
  free (msg_buf);
}


/*  MULTI_SUM: Global sums using multiple root nodes  */

void multi_sum (double *msg_buf, int buf_size)
{
  double msg_tmp[MAX_B];
  int    num_level, n_level, n_msg, n;
  int    span, span2, xchnge;
  int    send_loc, recv_loc;

  num_level = n_int(log((double) num_proc)/log(2.0));

  n = buf_size%num_proc;
  buf_size = buf_size/num_proc*num_proc;                         /* Set buffer to multiple of num_proc */
  if (n != 0)  global_sum (msg_buf+buf_size, n);                 /* Remnant by global_sum */
  buf_size = buf_size/num_proc;                                  /* Set buffer size */

  if ((1 << num_level) != num_proc)                              /* Use binary tree */
    {
      for (n = 0; n < num_proc; n++)
	global_sum (msg_buf+n*buf_size, buf_size);
    }
  else                                                           /* Use multi sum */
    {
      span = 1;  span2 = 2*span;                                 /* Shifts in processor and messages */
      for (n_level = 0; n_level < num_level; n_level++)          /* Gather and sum */
	{
	  xchnge = n_proc/span2*span2 + (n_proc + span)%span2;   /* Exchange processor id */
	  for (n_msg = 0; n_msg < num_proc/span2; n_msg++)
	    {
	      send_loc = ((xchnge + n_msg*span2)%num_proc)*buf_size;
	      recv_loc = ((n_proc + n_msg*span2)%num_proc)*buf_size;
	      vector_xchg (msg_buf+send_loc, msg_tmp, buf_size, xchnge, xchnge);
	      for (n = 0; n < buf_size; n++)
		msg_buf[recv_loc+n] += msg_tmp[n];
	    }
	  span *= 2;  span2 *= 2;                                /* Double span at each level */
	}

      for (n_level = 0; n_level < num_level; n_level++)          /* Scatter */
	{
	  span /= 2;  span2 /= 2;                                /* Halve span at each level */
	  xchnge = n_proc/span2*span2 + (n_proc + span)%span2;
	  for (n_msg = 0; n_msg < num_proc/span2; n_msg++)
	    {
	      send_loc = ((n_proc + n_msg*span2)%num_proc)*buf_size;
	      recv_loc = ((xchnge + n_msg*span2)%num_proc)*buf_size;
	      vector_xchg (msg_buf+send_loc, msg_buf+recv_loc, buf_size, xchnge, xchnge);
	    }
	}
    }
}
