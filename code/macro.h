/*  Macro header file  */
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


#define  n_int(x)       ( (x) >= (0) ? (int)((x)+0.5) : (int)((x)-0.5))
#define  s_int(x)       ( (x) >= (0) ? (int)(x) : (int)((x)-1.0))
#define  sgn(x)         ( (x) <  (0) ? (-1) : (1))
#define  max(a, b)      ( (a) >  (b) ? (a)  : (b))
#define  min(a, b)      ( (a) <  (b) ? (a)  : (b))
#define  box(a, b)      ( (a) - s_int((a)/(b))*(b))
#define  box2(a, b)     ( (a) >= (0.0) ? (a) - (int)((a)/(b))*(b): (a)-(int)((a)/(b)-1.0)*(b))
#define  n_image(a, b)  ( (a) - n_int((a)/(b))*(b))
#define  range(x, a, b) (((x) >= (a) && (x) < (b)) ? (1) : (0))
#define  cycle(x, a)  ( (x) >= 0 ? (x) % (a) : (x) % (a) + (a) )    // syfan 151101
