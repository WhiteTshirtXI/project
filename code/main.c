/*  Main routine for suspension code  */
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

int  num_proc, n_proc;
char *work_dir;
int main (int argc, char **argv)
{

//modified 20160916
//  char  *work_dir;
  
  init_procs (&argc, &argv);                          /* Initialization */
  num_proc = proc_num ();                             /* Number of procs */
  n_proc   = proc_id  ();                             /* Proc id */

  work_dir = argv[1];                                 /* Work directory */
 
  if(argc == 1)
    work_dir=".";
 
  if (n_proc == 0)
  {
    fprintf (stdout, "Running on %d processors\n", num_proc);
    fflush  (stdout);
  }
  
  double start = wclock();
  
  driver (work_dir);

  double finish = wclock();
  
  fini_procs ();                                      /* Wrap up */
  
  fprintf (stdout, "Elapsed time on proc %3d: %le (%le %le)\n", 
           n_proc, finish-start, start, finish);
  fflush  (stdout);

  return (0);
}
