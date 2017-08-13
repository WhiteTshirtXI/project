/*  Header File  */
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


/*  Define data types */

#define  Float  double

/*  Define hard limits */

#define  MAX_Y         2048          /* Max Y dimension */
#define  MAX_Z         1024          /* Max Z dimension */
#define  MAX_L         5000          /* Max # link cells */
#define  MAX_N          400          /* Max # neighbors */  // 200 modified
#define  MAX_BOND         6          /* Max # of bonds in blist */
#define  MAX_C        10000          /* Max # collisions */
#define  MAX_Cl        2000          /* Max cluster size for implicit update */
#define  MAX_W          100          /* Max # warnings */
#define  MAX_B         8192          /* Buffer size for msg passing */
#define  MAX_LENGTH    2000          /* Max chain length */
#define  MAX_COLL       500          /* Max # of colloids */
#define  MAX_DP         500          /* Max # of deformable particles */
#define  NTYPES           2          /* # of types of particles */

#define MAX_F            50

/*  Define parameters  */

#define  Num_Step         1          /* # steps for coord and vel update */
#define  Num_Hs3d_Step    1          /* # steps for hs3d update */
#define  Kappa          0.1          /* Stability Criteria for implicit update */
#define  Del_Sq         0.0          /* Max MSD before b_node update */
#define  Phi_Sed        1.0          /* Volume fraction at sedimentation front */
#define  Tol            1.0e-5       /* Tolerances */
#define  Elas_C         1.0          /* Elasticity of collisions */
#define  Num_Prop        44          /* # property vectors */
#define  Len             32          /* Integer word length */
#define  Pi            (acos(-1.0))
#define  MAX_EXT        0.9          /* Maximum chain extension */
#define  N_TYPES        2            /* # of particle types */

/* Define VSL RNG parameters */
#define BRNG    VSL_BRNG_MT19937
#define METHOD  VSL_METHOD_DGAUSSIAN_BOXMULLER2
#define SEED    1

/* Define output parameters */
#define WRITE_BUFFER     10          /* buffer for output */
#define BINS 40                      /* # of bins for the radial distribution function */
#define WRITE_BOND        0          /* write out bond connections */

/* Define pause command */
#define PAUSE printf("Press Enter key to continue..."); fgetc(stdin);

/*  Include other header files  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include "./lbe.h"
#include "./struct.h"
#include "./func_proto.h"
#include "./macro.h"
#include "mkl_vsl.h"
//#include "./mtwist.h"
//#include "./randistrs.h"
#include "./jacobi_eigenvalue.h"
