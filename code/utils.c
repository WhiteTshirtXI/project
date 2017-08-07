/*  UTILS:  Utility functions;
               FILE_NAME    (Constructs file name)
               WARNING      (Prints warning messages on processor)
               FATAL_ERR    (Fatal error handling)
               ERROR_CHK    (Fatal error checking)  */
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
static int  error_flag;


/*  FILE_NAME: Constructs file name labeled with task_number  */

void file_name (char *name, char *work_dir, int task_number)

{
  static char  name0[] = {'0','1','2','3','4','5','6','7','8','9'};
  int  index;
	
  if (work_dir != '\0')
  {
    for (index = strlen(name); index > 0; index--)
      name[index + strlen(work_dir)] = name[index-1];

    name[strlen(work_dir)] = '/';

    for (index = 0; index < strlen(work_dir); index++)
      name[index] = work_dir[index];
  }

  if(task_number >=0) 
  {
    index = strlen(name);

    name[index] = name0[task_number/100];

    task_number %= 100;

    name[index+1] = name0[task_number/10];

    task_number %= 10;

    name[index+2] = name0[task_number];

    name[index+3] = '\0';
  }
}


/*  WARNING: Prints warning messages on processor 0 */

void warning (char *warning_msg)
{
	static int n_warn=0;

	if (n_proc == 0)  fprintf (stdout, "WARNING: %s.\n", warning_msg);

	n_warn++;
	if (n_warn >= MAX_W)
	{
		fatal_err ("maximum number of warnings exceeded", -1);
		error_chk ();
	}
	fflush (stdout);
}


/*  FATAL_ERR: Fatal error handler  */

void fatal_err (char *error_msg, int flag)
{
 	error_flag = 1;
	if (flag < 0)
		fprintf (stdout, "FATAL ERROR: %s. Processor %d\n", error_msg, n_proc);
	else
		fprintf (stdout, "FATAL ERROR: %s %d. Processor %d\n", error_msg, flag, n_proc);
	fflush (stdout);
}


/*  ERROR_CHK Check for fatal errors */

void error_chk ()
{
	error_flag = global_max (error_flag);
	if (error_flag > 0)                          /* Fatal error: quit */
	{
		fprintf (stdout, "A FATAL ERROR HAS OCCURED: exiting proc %d\n", n_proc);
		fflush (stdout);
		exit(1);
	}      
}

void CheckVslError(int num)
{
    switch(num) {
        case VSL_ERROR_FEATURE_NOT_IMPLEMENTED: {
            printf("Error: this feature not implemented yet (code %d).\n",num);
            break;
        }
        case VSL_ERROR_UNKNOWN: {
            printf("Error: unknown error (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BADARGS: {
            printf("Error: bad arguments (code %d).\n",num);
            break;
        }
        case VSL_ERROR_MEM_FAILURE: {
            printf("Error: memory failure. Memory allocation problem maybe (code %d).\n",num);
            break;
        }
        case VSL_ERROR_NULL_PTR: {
            printf("Error: null pointer (code %d).\n",num);
            break;
        }
        case VSL_ERROR_INVALID_BRNG_INDEX: {
            printf("Error: invalid BRNG index (code %d).\n",num);
            break;
        }
        case VSL_ERROR_LEAPFROG_UNSUPPORTED: {
            printf("Error: leapfrog initialization is unsupported (code %d).\n",num);
            break;
        }
        case VSL_ERROR_SKIPAHEAD_UNSUPPORTED: {
            printf("Error: skipahead initialization is unsupported (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BRNGS_INCOMPATIBLE: {
            printf("Error: BRNGs are not compatible for the operation (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BAD_STREAM: {
            printf("Error: random stream is invalid (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BRNG_TABLE_FULL: {
            printf("Error: table of registered BRNGs is full (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BAD_STREAM_STATE_SIZE: {
            printf("Error: value in StreamStateSize field is bad (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BAD_WORD_SIZE: {
            printf("Error: value in WordSize field is bad (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BAD_NSEEDS: {
            printf("Error: value in NSeeds field is bad (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BAD_NBITS: {
            printf("Error: value in NBits field is bad (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BAD_UPDATE: {
            printf("Error: number of updated entries in buffer is invalid (code %d).\n",num);
            break;
        }
        case VSL_ERROR_NO_NUMBERS: {
            printf("Error: zero number of updated entries in buffer (code %d).\n",num);
            break;
        }
        case VSL_ERROR_INVALID_ABSTRACT_STREAM: {
            printf("Error: abstract random stream is invalid (code %d).\n",num);
            break;
        }
        case VSL_ERROR_FILE_CLOSE: {
            printf("Error: can`t close file (code %d).\n",num);
            break;
        }
        case VSL_ERROR_FILE_OPEN: {
            printf("Error: can`t open file (code %d).\n",num);
            break;
        }
        case VSL_ERROR_FILE_WRITE: {
            printf("Error: can`t write to file (code %d).\n",num);
            break;
        }
        case VSL_ERROR_FILE_READ: {
            printf("Error: can`t read from file (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BAD_FILE_FORMAT: {
            printf("Error: file format is unknown (code %d).\n",num);
            break;
        }
        case VSL_ERROR_UNSUPPORTED_FILE_VER: {
            printf("Error: unsupported file version (code %d).\n",num);
            break;
        }
    }

    if(num < 0) {
       exit(num);
    }
}


/* take the cross product c=(a)x(b) */
void product(double a[DIMS], double b[DIMS], double c[DIMS])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/* take the inner product c=(a)˙(b) */
double iproduct(double a[DIMS], double b[DIMS])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

