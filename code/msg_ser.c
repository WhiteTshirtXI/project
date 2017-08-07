/*  MSG_SER: Single processor "message-passing routines"; Contains
             VECTOR_COPY  (Copy a vector from send to recv)
             VECTOR_XCHG  (Exchange a vector with neighbors)
             BROAD_CAST   (Broadcast a vector)
             GLOBAL_SUM   (Global sum of a vector)
             GLOBAL_MAX   (Global maximum of a single integer)
             INIT_PROCS   (Initialize message passing)
             FINI_PROCS   (Finish with message passing)
             SYNC_PROCS   (Synchronize processors)
             NUM_PROC     (Get number of processors)
             N_PROC       (Get processor number)
             WCLOCK       (Wall clock timer)  */
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

/*  VECTOR_COPY: Copy a vector from send to recv */

void vector_copy (double *buf, int buf_siz, int n_proc, int send, int recv)
{
}


/*  VECTOR_XCHG: Receive from recv-Send to send */


void vector_xchg (double *send_buf, double *recv_buf, int buf_size, int send, int recv)
{       
	int  z;
	
	for (z = 0; z < buf_size; z++)
		recv_buf[z] = send_buf[z];
}


/*  GLOBAL_SUM: Sums contribution to a vector over all processors  */

void global_sum (double *vector, int buf_size)
{
}


/*  GLOBAL_MAX: Maximum of a single integer over all processors  */

int global_max (int g_max)
{
	return (g_max);
}


 /*  BROAD_CAST: Broadcast a vector to all processors  */

void broad_cast (double *vector, int buf_size, int n_proc)
{
}


/*  INIT_PROCS: Initialize message passing  */

void init_procs (int *argc, char ***argv)
{
}


/*  FINI_PROCS: Finish message passing  */

void fini_procs ()
{
}


/*  SYNC_PROCS: Synchronize processors  */

void sync_procs ()
{
}


/*  PROC_NUM: Get number of processors  */

int proc_num ()
{
	int  num_proc;

	num_proc = 1;
	return (num_proc);
}


/*  PROC_ID: Get processor number  */

int proc_id ()

{
	int  n_proc;

	n_proc = 0;
	return (n_proc);
}


/*  WCLOCK: Get time in secs  */

double wclock ()
{
#include <time.h>

	double  wclock;
	
	wclock  = (double) clock()/CLOCKS_PER_SEC;
	return (wclock);
}
