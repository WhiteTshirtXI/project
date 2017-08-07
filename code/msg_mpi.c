/*  MSG_MPI: MPI message-passing routines; Contains
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
#include  <mpi.h>


/*  VECTOR_COPY: Copy a vector from send to recv */

void vector_copy (double *buf, int buf_siz, int n_proc, int send, int recv)
{
	static MPI_Status   status;

	if (n_proc == send)  MPI_Send (buf, buf_siz, MPI_DOUBLE, recv, 0, MPI_COMM_WORLD);
	if (n_proc == recv)  MPI_Recv (buf, buf_siz, MPI_DOUBLE, send, 0, MPI_COMM_WORLD, &status);
}


/*  VECTOR_XCHG: Receive from recv-Send to send */

void vector_xchg (double *send_buf, double *recv_buf, int buf_siz, int send, int recv)
{
	static MPI_Request  request;
	static MPI_Status   status;

	MPI_Irecv(recv_buf, buf_siz, MPI_DOUBLE, recv, 0, MPI_COMM_WORLD, &request);
	MPI_Send (send_buf, buf_siz, MPI_DOUBLE, send, 0, MPI_COMM_WORLD);
	MPI_Wait (&request, &status);
}


/*  GLOBAL_SUM: Sums contribution to a vector over all processors  */

void global_sum (double *vector, int buf_size)

{
	double  work[MAX_B];
	int  i;

	MPI_Allreduce (vector, work, buf_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for (i = 0; i < buf_size; i++)
		vector[i] = work[i];
}


/*  GLOBAL_MAX: Maximum of an integer over all processors  */

int global_max (int g_max)
{
	int  work[1];
	
	MPI_Allreduce (&g_max, work, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	g_max = work[0];
	return (g_max);
}


/*  BROAD_CAST: Broadcast a vector to all processors  */

void broad_cast (double *vector, int buf_size, int n_proc)
{
  MPI_Bcast (vector, buf_size, MPI_DOUBLE, n_proc, MPI_COMM_WORLD);
}


/*  INIT_PROCS: Initialize message passing  */

void init_procs (int *argc, char ***argv)
{
	MPI_Init (argc, argv);
}


/*  FINI_PROCS: Finish message passing  */

void fini_procs ()
{
	MPI_Finalize ();
}


/*  SYNC_PROCS: Synchronize processors  */

void sync_procs ()
{
	MPI_Barrier (MPI_COMM_WORLD);
}


/*  PROC_NUM: Get number of processors  */

int proc_num ()
{
	int  num_proc;
	
	MPI_Comm_size (MPI_COMM_WORLD, &num_proc);
	return (num_proc);
}


/*  PROC_ID: Get processor number  */

int proc_id ()
{
	int  n_proc;
	
	MPI_Comm_rank (MPI_COMM_WORLD, &n_proc);
	return (n_proc);
}


/*  WCLOCK: Get time in secs  */

double wclock ()
{
	double wclock;
	
	wclock = MPI_Wtime();
	return (wclock);
}
