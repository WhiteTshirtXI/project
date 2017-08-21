/*  DRIVER:  Main driver routine for suspension code  */
/***********************************************************************
Susp3D: Lattice-Boltzmann simulation code for particle-fluid suspensions
Copyright (C) 2003 Tony Ladd

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is istributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 ***********************************************************************/

#include  "header.h"
#include  <time.h>

struct vector f_ext;
int    add_noise = FALSE;
double tau[2], tau_v[2], tau_g[2];
double lub_N, lub_T, lub_R;
double rho, ux, uy, uz;
double t_lbe;
int    max_x, max_y, max_z;
int    num_x, x_min, x_max;
int    num_sph, num_obj;
int    backflow_flag, wall_flag, shape_flag;
int    mass_flag, move_flag;
unsigned long  seed;
//mt_state twister;                        // random number generator (RNG) state
Float  std_xi[Num_Dir];                    // the amplitude of the fluctuation distributions
Float  mon_variance, mon_stdev;            // the amplitude of the monomer fluctuation
Float fluidflucvar[2], fluidflucstd[2];    // amplitude of fluid stress fluctuations 
Float std_xi[Num_Dir], std_xi_1[Num_Dir];  // amplitude of fluid stress fluctuations 

int initStep, growthStep; // Modification 20170712: Temporary design
double fictionalMass;
//
// Calls update and writes out averaged properties
// Declares node_map[XY][Z], velcs_df[XY][NDIR][Z]
//=================================================

void driver (char *work_dir)
{
	extern int    num_proc, n_proc;       // in our setup, num_proc=1 and n_proc=0 
	struct object *objects, *sph;
	struct vector f_sph, u_wall;
	struct sphere_param sphere_pm;
	struct sphere *spheres;
	struct monomer *monomers;
	struct face *faces;
	VSLStreamStatePtr rngstream;            // random number generator (RNG) state
	Float  ***velcs_df, **tmp_pp, *tmp_p;
	double y_prop[MAX_Y][Num_Prop]={{0}};
	double mass_fac, vel_fac, lub_cut, del_hy;
	double vol, min_rad=1e6, max_rad=0;
	double mon_variance;
	int    **node_map, *tmp_ip;
	int    i, j, k, nxy, q;
	int    size_xy, size_z;
	int    num_cycle, num_step, num_modes, n_cycle=0;
	int    n_sph, n_obj, n2_obj, n_mode, n, n_dir;
	int    task, num_slice, sum_slice, del_x;
	int    errcode;
	int    temp;
	char   task_file[200]={ "init/task_file.dat" };
	char   input_file[200]={ "init/input_file.dat" };
	char   alloc_file[200]={ "init/alloc_file.dat" };
	char   c_inp_file[200]={ "init/c_inp." };
	char   p_inp_file[200]={ "init/p_inp.dat" };
	char   properties[200]={ "data/properties.dat" };
	char   chk_p_file[200]={ "data/chk_p." };
	char   chk_f_file[200]={ "data/chk_f." };
	char   chk_c_file[200]={ "data/final.config" };
	char   p_out_file[200]={ "data/p_out." };
	char   m_out_file[200]={ "data/m_out." };
	char   modes_file[200]={ "data/modes." };
	FILE   *file_ptr=0;
	FILE   *stream;

	int mark_interval;
	int oscillation_period;
	int window;
	int flow_flag;
	int writeInterval_stress; // ctliao
  int retrieve_flow_step; // ctliao 20151110
	//double dqx, dqy, dqz;
	//int    num_qx, num_qy, num_qz;
  struct vector **forceDen, *pVector;

	// Set up file directory path 

	file_name (task_file, work_dir, -1);
	file_name (input_file, work_dir, -1);
	file_name (alloc_file, work_dir, -1);
	file_name (p_inp_file, work_dir, -1);
	file_name (properties, work_dir, -1);
	file_name (chk_c_file, work_dir, -1);

	// Get task

	file_ptr = fopen (task_file, "r");                               
	if (file_ptr == 0)  
		fatal_err("Could not open taskfile", -1);
	fscanf  (file_ptr, "%d", &task);
	fclose  (file_ptr);
	fprintf (stdout, "Begin driver: proc # %d, task %d\n", n_proc, task);
	fflush  (stdout);

	// Set up output file names

	file_name (c_inp_file, work_dir, task);
	file_name (chk_p_file, work_dir, task*num_proc+n_proc);
	file_name (chk_f_file, work_dir, task*num_proc+n_proc);
	file_name (p_out_file, work_dir, task*num_proc+n_proc);
	file_name (m_out_file, work_dir, task*num_proc+n_proc);
	file_name (modes_file, work_dir, task*num_proc+n_proc);

	// Read and print input data

	file_ptr = fopen (input_file, "r");
	if (file_ptr == 0)  
		fatal_err("Could not open input_file", -1);

	fscanf (file_ptr, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
	fscanf (file_ptr, "%d %d %d %d %d %d %le %d %d %le", &num_cycle, &num_step, &mark_interval, 
			&oscillation_period, &window, &num_modes, &t_lbe, &growthStep, &initStep, &fictionalMass);// Modification 20170712: Temporary design 
	fscanf (file_ptr, "%*s %*s %*s %*s");
	fscanf (file_ptr, "%d %d %d %d", &wall_flag, &flow_flag, &backflow_flag, &add_noise);
	fscanf (file_ptr, "%*s %*s %*s %*s");
	fscanf (file_ptr, "%le %le %le %le", &rho,  &ux, &uy, &uz);
	fscanf (file_ptr, "%*s %*s");
	fscanf (file_ptr, "%le %le",         &mass_fac,  &vel_fac);
	fscanf (file_ptr, "%*s %*s");
	fscanf (file_ptr, "%le %le",         &lub_cut,   &del_hy);
	fscanf (file_ptr, "%*s %*s %*s");
	fscanf (file_ptr, "%le %le %le",     &f_sph.x,   &f_sph.y,  &f_sph.z);
	fscanf (file_ptr, "%*s %*s %*s");
	fscanf (file_ptr, "%le %le %le",     &f_ext.x,   &f_ext.y,  &f_ext.z);
	fscanf (file_ptr, "%*s %*s %*s");
	fscanf (file_ptr, "%le %le %le",     &u_wall.x,  &u_wall.y, &u_wall.z);
	fscanf (file_ptr, "%*s %*s %*s");
	fscanf (file_ptr, "%le %le %le",     &tau[0],   &tau_v[0],    &tau_g[0]);
	fscanf (file_ptr, "%*s %*s %*s");
	fscanf (file_ptr, "%le %le %le",     &tau[1],   &tau_v[1],    &tau_g[1]);
	fscanf (file_ptr, "%*s %*s %*s");
	fscanf (file_ptr, "%le %le %le",     &lub_N,     &lub_T,    &lub_R);
	fscanf (file_ptr, "%*s ");
	fscanf (file_ptr, "%ld ",            &seed);
	fclose (file_ptr);

	srand( (unsigned)( time(NULL) ) );    // generate random seed
	seed = rand();

	if (n_proc == 0)                      // in our setup, n_proc = 0
	{
		file_ptr = fopen (properties, "w");
		if (file_ptr == 0)        fatal_err("Could not open properties", -1);

		fprintf (file_ptr, "Numbers (cycle, steps, modes_p, lbe_time, mark_interval, oscil_period, window): (%6d, %6d, %6d, %6lf, %d, %d, %d)\n",
				num_cycle, num_step, num_modes, t_lbe, mark_interval, oscillation_period, window);

		fprintf (file_ptr, "Flags (walls, backflow, add_noise): (%6d, %6d, %d)\n", wall_flag, 
				backflow_flag, add_noise);

		fprintf (file_ptr, "Rho, ux, uy, uz                (% .5e\t%.5e\t%.5e\t%.5e)\n",
				rho, ux, uy, uz);

		fprintf (file_ptr, "Mass fac, Vel_fac:             (% .5e\t% .5e)\n",  
				mass_fac,  vel_fac);

		fprintf (file_ptr, "Lub_cut   Del_hy:              (% .5e\t% .5e)\n",  
				lub_cut,   del_hy);

		fprintf (file_ptr, "Particle accelerations:        (% .5e\t% .5e\t% .5e)\n", 
				f_sph.x,  f_sph.y,  f_sph.z);

		fprintf (file_ptr, "External accelerations:        (% .5e\t% .5e\t% .5e)\n", 
				f_ext.x,  f_ext.y,  f_ext.z);

		fprintf (file_ptr, "Wall velocities:               (% .5e\t% .5e\t% .5e)\n", 
				u_wall.x, u_wall.y, u_wall.z);

		fprintf (file_ptr, "Relaxation times for fluid 1:  (% .5e\t% .5e\t% .5e)\n", 
				tau[0], tau_v[0], tau_g[0]);

		fprintf (file_ptr, "Relaxation times for fluid 2:  (% .5e\t% .5e\t% .5e)\n", 
				tau[1], tau_v[1], tau_g[1]);

		fprintf (file_ptr, "Lubrication ranges:            (% .5e\t% .5e\t% .5e)\n", 
				lub_N, lub_T, lub_R);

		fclose  (file_ptr);

	}

	/////////////////////////////////////////////////////////////////////////////////
	//                                                                             //
	// Flags:                                                                      //
	// wall_flag  0 = none;  1 = y;  2 = y & z;  3 = x & y & z;                    //
	// backflow_flag  0 = off;  1 = on;  2 = on + fluid at rest;                   //  
	//                3 = x only;  4 = y only;  5 = z only;                        //
	//                                                                             //
	// Forces: Input as accelerations                                              // 
	//                                                                             //
	// Y_prop:                                                                     //	
	// Fluid mass:             0                                                   //
	// Fluid momenta:	     1 = x,   2 = y,   3 = z                             //
	// Reynolds stress:	     4 = xx,  5 = yy,  6 = zz,  7 = yz,  8 = zx,  9 = xy //
	// Viscous stress:	    10 = xx, 11 = yy, 12 = zz, 13 = yz, 14 = zx, 15 = xy //
	// Particle mass:	    16                                                   //
	// Particle momenta:	    17 = ux, 18 = uy, 19 = uz, 20 = wx, 21 = wy, 22 = wz //
	// Particle stress:	    23 = xx, 24 = yy, 25 = zz, 26 = yz, 27 = zx, 28 = xy //
	// Collisional stress:    29 = xx, 30 = yy, 31 = zz, 32 = yz, 33 = zx, 34 = xy //
	// Particle-fluid stress: 35 = xx, 36 = yy, 37 = zz, 38 = yz, 49 = zx, 40 = xy //
	// Energy dissipation:    41 = t,  42 = r,  43 = c                             //
	//                                                                             //
	/////////////////////////////////////////////////////////////////////////////////

	/*  Read particle input files  ---------------------------------------------------------*/

	file_ptr = fopen (c_inp_file, "r");
	if (file_ptr == 0)        fatal_err("Could not open c_inp_file", -1); 

	fscanf(file_ptr, "%d %d %d %d", 
			&num_sph, &max_x, &max_y, &max_z);  /* # of spheres, maximum XYZ dimensions */

	if (max_x < 2)               fatal_err("max_x must be > 1", -1);
	if (max_y > MAX_Y)           fatal_err("array dimension MAX_Y must be > max_y", -1);
	if (max_z > MAX_Z)           fatal_err("array dimension MAX_Z must be > max_z", -1);
	if (MAX_B < max_z*Num_Dir_X) fatal_err("array dimension MAX_B must be > max_z*Num_Dir_X", -1);
	if (MAX_B < Num_Prop)        fatal_err("array dimension MAX_B must be > Num_Prop", -1);

	/* Memory for particles & walls */

	num_obj  = num_sph + wall_flag*2;

	objects  = (struct object *) calloc(num_obj, sizeof(struct object));   

	if (objects == 0)  
		fatal_err ("cannot allocate objects", -1);

	for (n_sph = 0; n_sph < num_sph; n_sph++)
	{
		sph = &objects[n_sph];

		/* read in particle position and orientation */

		fscanf (file_ptr, "%le %le %le %le %le %le", 
				&sph->r.x, &sph->r.y, &sph->r.z, &sph->e.x, &sph->e.y, &sph->e.z);

		/* read in paricle linear and angular velocity */

		fscanf (file_ptr, "%le %le %le %le %le %le", 
				&sph->u.x, &sph->u.y, &sph->u.z, &sph->w.x, &sph->w.y, &sph->w.z);

		//////////////////////////////////////////////////////////////
		//                                                          //
		//  particle shape: 0=sphere, 1=ellipsoid, 2=open cylinder, // 
		//                  3=closed cylinder, 4=capped cylinder    // 
		//  particle mass:  0=infinity, 1=finite mass               //
		//                                                          //
		//  particle move:  0=fixed particle, 1=standard move,      // 
		//                  2=particle as part of wall              //
		//                                                          //
		//////////////////////////////////////////////////////////////

		/*  particle density and radii */

		fscanf (file_ptr, "%d %d %d %le %le %le", 
				&sph->shape_flag, &sph->mass_flag, 
				&sph->move_flag, &sph->rho, &sph->r_a, &sph->r_b);

		min_rad = min(min_rad, min(sph->r_a,sph->r_b));
		max_rad = max(max_rad, max(sph->r_a,sph->r_b));
	}
	fclose (file_ptr);
	file_ptr = 0;

	if (min_rad < sqrt(1.25))
		warning ("smallest particle dimension should be at least sqrt(1.25)");
	if (min(max_x, min(max_y,max_z)) <= max_rad*4 && num_sph > 1)
		warning ("smallest box dimension should be at least 4 times the maximum radius");
	if (min(max_x, min(max_y,max_z)) <= max_rad*2+1 && num_sph == 1)
		warning ("smallest box dimension should be at least 2 times the maximum radius");

	for(n_sph = 0; n_sph < num_sph; n_sph ++)
	{
		printf("Particle %d : r=(%g %g %g), u=(%g %g %g), rho=%le, R=%le\n", n_sph, 
				objects[n_sph].r.x, objects[n_sph].r.y, objects[n_sph].r.z, 
				objects[n_sph].u.x, objects[n_sph].u.y, objects[n_sph].u.z, 
				objects[n_sph].rho, objects[n_sph].r_a);
	}

	// Read point force polymer input file

	file_ptr = fopen (p_inp_file, "r");
	if (file_ptr == 0)    fatal_err("Could not open p_inp.dat", -1);

	//////////////////////////////////////////////////////////////////
	//                                                              //
	// # of spheres; maximum XYZ dimensions;                        //
	// spring = 0 (FENE), 1(WLC);                                   //
	// ev_type = 0 (HS), 1(WCA), 2(gaussian);                       //
	// verlet_type=1 (2nd order), 0 (first order);                  //
	// Spring params:                                               //
	// H_FENE, Q_FENE, kuhn length, # of kuhn segments per spring,  // 
	// bending strength;                                            //
	// evcutoff = exc. vol. cutoff parameter;                       //
	// fric = friction coef. of a point force;                      //
	// dt = timestep;                                               //
	//                                                              //
	//////////////////////////////////////////////////////////////////

	fscanf(file_ptr, "%*s %*s %*s %*s %*s %*s %*s %*s ");
	fscanf(file_ptr, "%d %d %d %d %d %d %d", 
			&sphere_pm.Ntype[0], &sphere_pm.Ntype[1], 
			&sphere_pm.nlevel[0], &sphere_pm.nlevel[1],   &max_x, &max_y, &max_z);  
	fscanf(file_ptr, "%*s %*s %*s %*s %*s");
	fscanf(file_ptr, "%d %d %d %d", &sphere_pm.spring, 
			&sphere_pm.ev_type, &sphere_pm.verlet, &sphere_pm.initconfig);
	fscanf(file_ptr, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");  
	fscanf(file_ptr, "%le %le %le %le %le %le %le %le %le %le %le %le %le",
			&sphere_pm.H_fene[0], &sphere_pm.H_fene[1], 
			&sphere_pm.Q_fene[0], &sphere_pm.Q_fene[1], 
			&sphere_pm.sigma_k, &sphere_pm.nks, 
			&sphere_pm.k_bend[0], &sphere_pm.k_bend[1], 
			&sphere_pm.k_V[0], &sphere_pm.k_V[1], 
			&sphere_pm.k_A[0], &sphere_pm.k_A[1], &sphere_pm.ka_local);
	fscanf(file_ptr, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
	fscanf(file_ptr, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le", 
      &sphere_pm.evcutoff,&sphere_pm.fric,&sphere_pm.dt,&sphere_pm.monmass,&sphere_pm.kT, 
      &sphere_pm.epsilon_LJ,&sphere_pm.cutoff_LJ,&sphere_pm.epsilon_LJ_final, 
      &sphere_pm.cutoff_LJ_final,&sphere_pm.eq_LJ,&sphere_pm.epsilon_WCA,&sphere_pm.eq_WCA, 
      &sphere_pm.epsilon_overlap,&sphere_pm.eq_overlap);
  fscanf(file_ptr, "%*s %*s %*s %*s %*s %*s %*s");
  fscanf(file_ptr, "%d %d %le %le %le %le", &sphere_pm.attrac_type,&sphere_pm.attrac_type_1,
      &sphere_pm.eq_morse,&sphere_pm.depth_morse,&sphere_pm.width_morse,
      &sphere_pm.cutoff_morse);
	fscanf(file_ptr, "%*s %*s");
	fscanf(file_ptr, "%le %le %le", 
			&sphere_pm.f_ext[0], &sphere_pm.f_ext[1], &sphere_pm.f_ext[2]);
	fscanf(file_ptr, "%*s %*s %*s %*s %*s");
	fscanf(file_ptr, "%d %d %d %d", &sphere_pm.write_time, &sphere_pm.write_config, 
			&sphere_pm.write_fluid, &writeInterval_stress);
  fscanf(file_ptr,"%*s %*s");
  fscanf(file_ptr,"%d",&retrieve_flow_step);

	if (max_x < 2)               fatal_err("max_x must be > 1", -1);
	if (max_y > MAX_Y)           fatal_err("array dimension MAX_Y must be > max_y", -1);
	if (max_z > MAX_Z)           fatal_err("array dimension MAX_Z must be > max_z", -1);
	if (MAX_B < max_z*Num_Dir_X) fatal_err("array dimension MAX_B must be > max_z*Num_Dir_X", -1);
	if (MAX_B < Num_Prop)        fatal_err("array dimension MAX_B must be > Num_Prop", -1);
	fclose(file_ptr);

	// Allocation of slices to each processor ??

	file_ptr = fopen(alloc_file, "r");  
	if(file_ptr != 0) 
  {
		sum_slice = 0;
		for(n = 0; n < num_proc; n++) {
			fscanf (file_ptr, "%d", &num_slice);
			sum_slice += num_slice;      
			if(n == n_proc) {
				num_x = num_slice;
				x_min = sum_slice - num_x;
				x_max = sum_slice;
			}
		}
		fprintf(stdout, "Allocation to proc %3d: %4d %4d %4d\n", n_proc, num_x, x_min, x_max);       
    fflush (stdout);    
		if(sum_slice != max_x)    fatal_err("incorrect allocation of slices", -1);
	}
	else  // run this block.
	{
		num_x = max_x / num_proc;                      /* # of slices per processor */
		del_x = max_x % num_proc;                      /* # of extra slices */
		x_min = n_proc * num_x + min(n_proc,del_x);    /* Processor boundaries */
		x_max = x_min + num_x;
		if(del_x > n_proc)    x_max++;
		if(del_x != 0)    warning ("number of yz-slices not divisible by num_proc");
	}
	if(num_x == 1)    fatal_err("only one yz-slice", -1);

	//  Get space and set pointers for population densities and node map

	size_xy = (num_x+2) * (max_y+2);  /* Allow for end slices */
	size_z  = max_z+2;
	velcs_df = (Float ***) calloc (size_xy, sizeof(*velcs_df));
	if(velcs_df == 0)  
	  fatal_err("cannot allocate population densities", -1);			   
	tmp_pp = (Float **) calloc (size_xy*Num_Dir, sizeof(*tmp_pp));
	if(tmp_pp == 0)  
	  fatal_err("cannot allocate population densities", -1);
	for(n = 0; n < size_xy; n++) {   /* Assign pointers to pointers */
	  velcs_df[n] = tmp_pp;
		tmp_pp += Num_Dir;
	}
	tmp_p = (Float *) calloc (size_xy*Num_Dir*size_z, sizeof(*tmp_p));
	if(tmp_p == 0)  
	  fatal_err("cannot allocate population densities", -1);
	for(n = 0; n < size_xy; n++) {
		for(n_dir = 0; n_dir < Num_Dir; n_dir++) {
			velcs_df[n][n_dir] = tmp_p;
			tmp_p += size_z;
		}
	}
  /* Note 20170820 
  -----------------
  Assign initial values to the distribution function. 
  The initial fluid velocity is choosen to be 0, and the initial density is choosen to be 
  36. To freely assign a value for the fluid density variable fac[19] in lbe.h should be 
  modified and the moving boundary condition in bnodes.c should be tooken care.
  */
  double density = 36.0;
  for(i=1; i <= max_x; i++) {
    for(j=1; j <= max_y; j++) {
      nxy = i*(max_y+2) + j;
      for(k=1; k <= max_z; k++) {
        for(q=0; q < 19; q++) {
          velcs_df[nxy][q][k] = weight[q]*density; 
        }
      }
    }
  }

	node_map = (int **) calloc (size_xy, sizeof(*node_map));
	tmp_ip   = (int *)  calloc (size_xy*size_z, sizeof(*tmp_ip));
	for(n = 0; n < size_xy; n++) {
		node_map[n] = tmp_ip;
		tmp_ip += size_z;
	}

  forceDen = (struct vector **) calloc(size_xy, sizeof(struct vector *)); // Modification 20170419
  if(forceDen == 0) { fatal_err("cannot allocate forceDen", -1); }
  pVector = (struct vector *) calloc(size_xy*size_z, sizeof(struct vector));
  if(pVector == 0) { fatal_err("cannot allocate pVector of forceDen",-1); }
  for(int i=0; i < size_xy; i++) {
    forceDen[i] = pVector;
    pVector += size_z;
  }
  pVector = NULL;

	error_chk();    /* Exit if fatal error in startup */

	/*//////////////////////////////////////////////////////////////////////////////////////*/

	for (j=0 ; j<NTYPES ; j++) 
	{
		if(sphere_pm.nlevel[j] >= 0) 
		{
			sphere_pm.N_per_sphere[j] = 12;
			sphere_pm.face_per_sphere[j] = 20;
			temp = 30;

			for (i=0 ; i<sphere_pm.nlevel[j] ; i++) 
			{
				sphere_pm.N_per_sphere[j] += temp;
				sphere_pm.face_per_sphere[j] *= 4;
				temp *= 4;
			}
		}
		else if(sphere_pm.nlevel[j] == -1) 
		{
			sphere_pm.N_per_sphere[j] = 162;
			sphere_pm.face_per_sphere[j] = 320;
		}
	}

	sphere_pm.Nsphere   = sphere_pm.Ntype[0] + sphere_pm.Ntype[1];

	sphere_pm.num_beads = sphere_pm.Ntype[0]*sphere_pm.N_per_sphere[0] + 
		sphere_pm.Ntype[1]*sphere_pm.N_per_sphere[1];

	sphere_pm.nfaces    = sphere_pm.Ntype[0]*sphere_pm.face_per_sphere[0] + 
		sphere_pm.Ntype[1]*sphere_pm.face_per_sphere[1];

	sphere_pm.numsteps  = num_step*num_cycle;

	sphere_pm.MD_steps  = (int)((t_lbe+1e-6)/sphere_pm.dt);  /* # of MD steps 
																															per LB step */
	/* variance of the monomer fluctuation force */
	mon_variance = 2.0*sphere_pm.kT*Rho_Fl*sphere_pm.fric*(1.0/sphere_pm.dt);
	mon_stdev = sqrt(mon_variance);

	/*sphere_pm.sigma_k = 0.106/(0.077/(sphere_pm.fric / (Pi*(2.0*tau-1.0))));*/

	sphere_pm.Ss = sphere_pm.sigma_k*sqrt(sphere_pm.nks/6.0);

	sphere_pm.relax_time = 10000;

	/* Ladd JFM 1994 */
	fluidflucvar[0] = Rho_Fl*sphere_pm.kT / 3.0 *(1.0-(1.0-1.0/tau[0])*(1.0-1.0/tau[0])) ;
	fluidflucstd[0] = sqrt(fluidflucvar[0]);
	fluidflucvar[1] = Rho_Fl*sphere_pm.kT / 3.0 *(1.0-(1.0-1.0/tau[1])*(1.0-1.0/tau[1])) ;
	fluidflucstd[1] = sqrt(fluidflucvar[1]);

	for(q=4; q<7; q++)
		std_xi[q]= fluidflucstd[0]*sqrt(4.0/3.0);
	for(q=7; q<10; q++)
		std_xi[q]= fluidflucstd[0];

	for(q=4; q<7; q++)
		std_xi_1[q]= fluidflucstd[1]*sqrt(4.0/3.0);
	for(q=7; q<10; q++)
		std_xi_1[q]= fluidflucstd[1];

	/* mts_seed32new( &twister, seed);  */

	/* Initialize the fluctuation generator */
	errcode=vslNewStream( &rngstream, BRNG,  seed);
	CheckVslError(errcode);

	Float omega = 1.0/tau[0];
	Float ghost_omega = 1.0/tau_g[0];

	//  Initialize spheres

	if(sphere_pm.Nsphere > 0) 
  {
		spheres = (struct sphere*)malloc(sphere_pm.Nsphere*sizeof(struct sphere));
		if(spheres == NULL) {
			printf("Memory allocation for spheres failed\n");
			exit(1);
		}
		monomers = (struct monomer*)malloc(sphere_pm.num_beads*sizeof(struct monomer));
		if(monomers == NULL) {
			printf("Memory allocation for monomers failed\n");
			exit(2);
		}
		faces = (struct face*)malloc(sphere_pm.nfaces*sizeof(struct face));
		if(faces == NULL) {
			printf("Memory allocation for monomers failed\n");
			exit(2);
		}
		sphere_init (&sphere_pm, spheres, monomers, faces, work_dir); 
    read_particle_para (&sphere_pm, monomers, faces, work_dir);  // Modification 0422
    calculate_particle_para (&sphere_pm, monomers, faces, work_dir);
		bnodes_dp (node_map, &sphere_pm, spheres, monomers);
	}

	// Initialize objects

	f_ext.x *= Rho_Fl;  
	f_ext.y *= Rho_Fl;  
	f_ext.z *= Rho_Fl;
	vol = max_x * max_y * max_z;

	/* Spheroidal Particles */

	for (n_sph = 0; n_sph < num_sph; n_sph++)           
	{
		objects[n_sph].lub_cut   = lub_cut;
		objects[n_sph].del_hy    = del_hy;

		spheroid_init (&objects[n_sph], mass_fac, vel_fac);  /* create spheres at given 
																														positions and velocity */   
		if (objects[n_sph].mass_flag)
		{
			objects[n_sph].f_ext.x = f_sph.x*objects[n_sph].mass;
			objects[n_sph].f_ext.y = f_sph.y*objects[n_sph].mass;
			objects[n_sph].f_ext.z = f_sph.z*objects[n_sph].mass;

			if (backflow_flag)
			{
				if (backflow_flag != 4 && backflow_flag != 5) 
					objects[n_sph].f_ext.x *= 1.0 - 1.0/(objects[n_sph].rho*mass_fac);
				if (backflow_flag != 3 && backflow_flag != 5) 
					objects[n_sph].f_ext.y *= 1.0 - 1.0/(objects[n_sph].rho*mass_fac);
				if (backflow_flag != 3 && backflow_flag != 4) 
					objects[n_sph].f_ext.z *= 1.0 - 1.0/(objects[n_sph].rho*mass_fac);
			}

			if (backflow_flag)
			{
				f_ext.x -= objects[n_sph].f_ext.x/vol;
				f_ext.y -= objects[n_sph].f_ext.y/vol;
				f_ext.z -= objects[n_sph].f_ext.z/vol;
			}

		}
	}

	sphere_pm.f_ext[0] *= Rho_Fl; 
	sphere_pm.f_ext[1] *= Rho_Fl; 
	sphere_pm.f_ext[2] *=Rho_Fl;

	/* monomers */

	for (n_sph = 0; n_sph < sphere_pm.num_beads; n_sph++)          
	{
		monomers[n_sph].f_ext[0] = sphere_pm.f_ext[0]*monomers[n_sph].rho;
		monomers[n_sph].f_ext[1] = sphere_pm.f_ext[1]*monomers[n_sph].rho;
		monomers[n_sph].f_ext[2] = sphere_pm.f_ext[2]*monomers[n_sph].rho;

		if (backflow_flag) 
		{
			f_ext.x -= monomers[n_sph].f_ext[0]/vol;
			f_ext.y -= monomers[n_sph].f_ext[1]/vol;
			f_ext.z -= monomers[n_sph].f_ext[2]/vol;
		}
	}

	if (backflow_flag)  
		warning ("backflow forces were added to external force density");

	/* initialization of RNG for each proc > 0 */
	if(n_proc > 0)
		/* mts_seed32new( &twister, seed+n_proc*161287); */
		vslNewStream( &rngstream, BRNG,  seed%n_proc*(seed/n_proc));

	/*  Set up container walls; only certain flows 
			are possible for each wall_flag setting (see annotations)  */

	for (n2_obj = num_sph; n2_obj < num_obj; n2_obj++)     /* Lubrication cutoff */
		objects[n2_obj].lub_cut = lub_cut;

	n_obj = num_sph;

	if (wall_flag == 1) /* parallel slit channel */
	{
		wall_init (&objects[n_obj], 0);           /* Y- Wall */  // ctliao 20151110
		objects[n_obj].u_init.x = -0.5*u_wall.x;  /* Planar shear */  // 20161227
		objects[n_obj].u_init.y =      u_wall.y;  /* Fluidization */
		objects[n_obj].u_init.z =      u_wall.z;  /* Channel flow */
    objects[n_obj].u.x = 0.0; 
    objects[n_obj].u.y = 0.0;
    objects[n_obj].u.z = 0.0; 
		n_obj++;

		wall_init (&objects[n_obj], 1);                    /* Y+ Wall */
		objects[n_obj].u_init.x =  0.5*u_wall.x;  //20161227
		objects[n_obj].u_init.y =      u_wall.y;
		objects[n_obj].u_init.z =      u_wall.z;
    objects[n_obj].u.x = 0.0;
    objects[n_obj].u.y = 0.0;
    objects[n_obj].u.z = 0.0;
		n_obj++;

		if  (del_hy != 0.0)
		{
			warning ("particles must not be within del_hy of the y-wall");
			objects[n_obj-2].r.y += del_hy;    // Adjust wall location
			objects[n_obj-1].r.y -= del_hy;    // for lubrication
		}

		if  (u_wall.y != 0.0)
		{
			warning ("particles must not be within 1 lattice unit of the y-wall");
			objects[n_obj-2].r.y ++;    // Adjust wall location for lubrication 
			objects[n_obj-1].r.y --;
		}
	}

	if (wall_flag == 2)  /* rectangular channel */
	{
		wall_init (&objects[n_obj], 0);                    /* Y- Wall */
		objects[n_obj].u.x =      u_wall.x;                /* Channel flow */
		objects[n_obj].u.y =  0.5*u_wall.y;                /* Elongational flow */
		n_obj++;

		wall_init (&objects[n_obj], 1);                    /* Y+ Wall */
		objects[n_obj].u.x =      u_wall.x;
		objects[n_obj].u.y = -0.5*u_wall.y;
		n_obj++;

		wall_init (&objects[n_obj], 2);                    /* Z- Wall */
		objects[n_obj].u.x =      u_wall.x;
		objects[n_obj].u.z = -0.5*u_wall.y;
		n_obj++;

		wall_init (&objects[n_obj], 3);                    /* Z+ Wall */
		objects[n_obj].u.x =      u_wall.x;
		objects[n_obj].u.z =  0.5*u_wall.y;
		n_obj++;

		if (del_hy != 0.0)
		{
			warning ("particles must not be within del_hy of the y-wall");
			warning ("particles must not be within del_hy of the z-wall");
			objects[n_obj-4].r.y += del_hy;
			objects[n_obj-3].r.y -= del_hy;
			objects[n_obj-2].r.z += del_hy;
			objects[n_obj-1].r.z -= del_hy;
		}
		if (u_wall.y != 0.0)
		{
			if (max_y != max_z)  fatal_err ("max_y & max_z must be equal for this flow", -1);
			warning ("z-wall velocity has been set to -y-wall velocity");
			warning ("particles must not be within 1 lattice unit of the y-wall");
			warning ("particles must not be within 1 lattice unit of the z-wall");
			objects[n_obj-4].r.y ++;
			objects[n_obj-3].r.y --;
			objects[n_obj-2].r.z ++;
			objects[n_obj-1].r.z --;
		}
	}

	if (wall_flag == 3)   /* boxed geometry */
	{
		wall_init (&objects[n_obj], 0);                    /* Y- Wall */
		objects[n_obj].u.y =      u_wall.y;                /* Fluidization */
		n_obj++;
		wall_init (&objects[n_obj], 1);                    /* Y+ Wall */
		objects[n_obj].u.y =      u_wall.y;
		n_obj++;
		wall_init (&objects[n_obj], 2);                    /* Z- Wall */
		n_obj++;
		wall_init (&objects[n_obj], 3);                    /* Z+ Wall */
		n_obj++;
		wall_init (&objects[n_obj], 4);                    /* X- Wall */
		n_obj++;
		wall_init (&objects[n_obj], 5);                    /* X+ Wall */
		n_obj++;

		if (del_hy != 0.0)
		{
			warning ("particles must not be within del_hy of the x-wall");
			warning ("particles must not be within del_hy of the y-wall");
			warning ("particles must not be within del_hy of the z-wall");
			objects[n_obj-6].r.y += del_hy;
			objects[n_obj-5].r.y -= del_hy;
			objects[n_obj-4].r.z += del_hy;
			objects[n_obj-3].r.z -= del_hy;
			objects[n_obj-2].r.z += del_hy;
			objects[n_obj-1].r.z -= del_hy;
		}
		if  (u_wall.x != 0.0)  warning ("x-wall velocity has been set to 0");
		if  (u_wall.y != 0.0)
		{
			warning ("particles must not be within 1 lattice unit of the y-wall");
			objects[n_obj-6].r.y ++;    // Adjust wall location for lubrication 
			objects[n_obj-5].r.y --;
		}
		if  (u_wall.z != 0.0)  warning ("z-wall velocity has been set to 0");
	}

	/* Read in external particle and fluid forces -------------------------------------------*/

	if (n_proc == 0)
	{ 
		file_ptr = fopen (properties, "a");   
		if (file_ptr == 0)        fatal_err("Could not open properties", -1); 

		fprintf (file_ptr, "\nExternal particle forces:\n");

		for (n_sph = 0; n_sph < num_sph; n_sph++)
			fprintf (file_ptr, "%6d % .5e % .5e % .5e\n", n_sph, 
					objects[n_sph].f_ext.x, objects[n_sph].f_ext.y, objects[n_sph].f_ext.z);

		fprintf (file_ptr, "\nExternal fluid forces:\n");
		fprintf (file_ptr, "% .5e % .5e % .5e\n", 
				f_ext.x*vol, f_ext.y*vol, f_ext.z*vol);

		fclose  (file_ptr);
	}

	if(n_proc == 0) 
	{
		file_ptr = fopen(properties, "a");

		fprintf(file_ptr, "num_walls = %d \n", n_obj-n_sph);

		fprintf(file_ptr, "kinematic viscosity = (%le; %le)\n", 
				(tau[0]-0.5)/3.0, (tau[1]-0.5)/3.0);

		fprintf(file_ptr,"\ntotal steps = %d \nMD timestep = %le, %d MD steps per LBE step\nseed =% d\n", 
				sphere_pm.numsteps, sphere_pm.dt, sphere_pm.MD_steps, seed);

		fprintf(file_ptr,"%d spheres\n%d num_beads\n%d num_faces\n verlet type%d\n spring type%d\n ev type %d\n",
				sphere_pm.Nsphere, sphere_pm.num_beads, sphere_pm.nfaces, 
				sphere_pm.verlet, sphere_pm.spring, sphere_pm.ev_type);

		fprintf(file_ptr, "H_fene=(%le %le), Q_fene=(%le %le)\n", 
				sphere_pm.H_fene[0], sphere_pm.H_fene[1], 
				sphere_pm.Q_fene[0], sphere_pm.Q_fene[1]);

		fprintf(file_ptr, "k_bend=(%le %le), k_V=(%le %le), k_A=(%le %le)\n", 
				sphere_pm.k_bend[0], sphere_pm.k_bend[1], 
				sphere_pm.k_V[0], sphere_pm.k_V[1], 
				sphere_pm.k_A[0], sphere_pm.k_A[1]);

		fprintf(file_ptr, "evcutoff=%le, fric=%le, mass=%le, dt=%le\n", 
				sphere_pm.evcutoff, sphere_pm.fric, sphere_pm.monmass, sphere_pm.dt);

		fprintf(file_ptr, "kT = %le monomer fluc stdev = %le, fluid fluc var=%le %le\n", 
				sphere_pm.kT, mon_stdev, fluidflucvar[0], fluidflucvar[1]);

		fprintf(file_ptr, "monomer hydrodynamic radius = %le, Ss=%le\n", 
				sphere_pm.fric/(Pi*(2.0*tau[0]-1.0)), sphere_pm.Ss);

		fprintf(file_ptr, "sphere relaxation time = %d\n", sphere_pm.relax_time);

		fprintf(file_ptr, "omega = %le, ghost_omega = %le\n", omega, ghost_omega);

		fprintf(file_ptr, "tau = (%le %le), tau_v= (%le %le), tau_g = (%le %le)\n", 
				tau[0], tau[1], tau_v[0], tau_v[1], tau_g[0], tau_g[1]);

		fprintf(file_ptr, "extern monomer force = %le %le %le\n", 
				sphere_pm.f_ext[0], sphere_pm.f_ext[1], sphere_pm.f_ext[2]);

		fprintf(file_ptr, "particle V0= (%le %le) A0=(%le %le)\n", 
				sphere_pm.V0[0], sphere_pm.V0[1], sphere_pm.A0[0], sphere_pm.A0[1]);

		fprintf(file_ptr, "output props every %d steps, config every %d steps, fluid every %d steps\n", 
				sphere_pm.write_time, sphere_pm.write_config, sphere_pm.write_fluid);

		fclose(file_ptr);
	}

	if (n_proc == 0)
	{
		file_ptr = fopen (properties, "a");
		if (file_ptr == 0)  
			fatal_err("Could not open properties", -1);
		fclose  (file_ptr);
	}

	// Main Loop 
	while (n_cycle < num_cycle)                           
	{
		//--------------------
		//  Run the LBE update 
		//--------------------
		n_cycle = update (num_cycle, num_step, num_modes, flow_flag, mark_interval, 
	      oscillation_period, window, backflow_flag, rngstream, node_map, velcs_df, y_prop, 
        modes_file, chk_p_file, chk_f_file, chk_c_file, f_ext, objects, spheres, monomers, 
				faces, &sphere_pm, work_dir, n_cycle, writeInterval_stress, retrieve_flow_step, 
        forceDen);
		//-----------------------  
		//  Normalize output data
		//-----------------------    
		for (n_obj = 0; n_obj < num_obj; n_obj++)
		{
			objects[n_obj].p.x /= (double) num_step;          /*  Average forces  */
			objects[n_obj].p.y /= (double) num_step;
			objects[n_obj].p.z /= (double) num_step;
			objects[n_obj].l.x /= (double) num_step;          /*  Average torques  */
			objects[n_obj].l.y /= (double) num_step;
			objects[n_obj].l.z /= (double) num_step;
		}

		for (n = 1; n <= max_y; n++)
			for (n_mode = 0; n_mode < Num_Prop; n_mode++)
				y_prop[n][n_mode] /= max_x*max_z*num_step;

		//-------------------
		//  Write output data  
		//-------------------    
		if ((n_cycle-1)%num_proc == n_proc)
		{
			fprintf (stdout, "Begin data output  %5d on %3d\n", n_cycle, n_proc);
			fflush (stdout);

			file_ptr = fopen (p_out_file, "a");
			if (file_ptr == 0)  fatal_err("Could not open p_out_file", -1);
			fflush (file_ptr);

			fprintf (file_ptr, "Cycle %6d: Coordinates\n", n_cycle);

			for (n_sph = 0; n_sph < num_sph; n_sph++)
				fprintf (file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", 
						n_sph, objects[n_sph].r.x, objects[n_sph].r.y, 
						objects[n_sph].r.z, objects[n_sph].e.x, objects[n_sph].e.y, 
						objects[n_sph].e.z);

			fprintf (file_ptr, "      %6d: Velocities\n", n_cycle);

			for (n_sph = 0; n_sph < num_sph; n_sph++)
				fprintf (file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", 
						n_sph, objects[n_sph].u.x, objects[n_sph].u.y, 
						objects[n_sph].u.z, objects[n_sph].w.x, objects[n_sph].w.y, 
						objects[n_sph].w.z);

			fprintf (file_ptr, "      %6d: Forces\n", n_cycle);

			for (n_sph = 0; n_sph < num_sph; n_sph++)
				fprintf (file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", 
						n_sph, objects[n_sph].p.x, objects[n_sph].p.y, 
						objects[n_sph].p.z, objects[n_sph].l.x, objects[n_sph].l.y, 
						objects[n_sph].l.z);

			fprintf (file_ptr, "      %6d: Average Wall Forces\n", n_cycle);

			for (n_obj = num_sph; n_obj < num_obj; n_obj++)
				fprintf (file_ptr, " % .5e % .5e % .5e\n",
						objects[n_obj].p.x, objects[n_obj].p.y, objects[n_obj].p.z);

			fclose (file_ptr);

			file_ptr = fopen(m_out_file, "a");
			if (file_ptr == 0)  fatal_err("Could not open m_out_file", -1);
			fflush (file_ptr);

			fprintf (file_ptr, "Cycle %6d: Fluid Mass/Momenta\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e % .5e\n", n,
						y_prop[n][ 0], y_prop[n][ 1], y_prop[n][ 2], y_prop[n][ 3]);

			fprintf (file_ptr, "      %6d: Reynolds Pressure\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", n,
						y_prop[n][ 4], y_prop[n][ 5], y_prop[n][ 6], y_prop[n][ 7], 
						y_prop[n][ 8], y_prop[n][ 9]);

			fprintf (file_ptr, "      %6d: Viscous Pressure\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", n,
						y_prop[n][10], y_prop[n][11], y_prop[n][12], y_prop[n][13], 
						y_prop[n][14], y_prop[n][15]);

			fprintf (file_ptr, "      %6d: Particle Mass/Momenta\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e % .5e\n",
						n, y_prop[n][16], y_prop[n][17], y_prop[n][18], 
						y_prop[n][19], y_prop[n][20], y_prop[n][21], y_prop[n][22]);

			fprintf (file_ptr, "      %6d: Particle Pressure\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", 
						n, y_prop[n][23], y_prop[n][24], y_prop[n][25], 
						y_prop[n][26], y_prop[n][27], y_prop[n][28]);

			fprintf (file_ptr, "      %6d: Collisional Pressure\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", n,
						y_prop[n][29], y_prop[n][30], y_prop[n][31], y_prop[n][32], 
						y_prop[n][33], y_prop[n][34]);

			fprintf (file_ptr, "      %6d: Particle-Fluid Pressure\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e % .5e % .5e % .5e\n", n,
						y_prop[n][35], y_prop[n][36], y_prop[n][37], y_prop[n][38], 
						y_prop[n][39], y_prop[n][40]);

			fprintf (file_ptr, "      %6d: Energy dissipation\n", n_cycle);
			for (n = 1; n <=  max_y; n++)
				fprintf(file_ptr, "%6d % .5e % .5e % .5e\n", n,
						y_prop[n][41], y_prop[n][42], y_prop[n][43]);
			fclose (file_ptr);

			fprintf (stdout, "End data output    %5d on %3d\n", n_cycle, n_proc);
			fflush (stdout);
		}

		//------------------
		// Zero output data  
		//------------------	   
		for (n = 1; n <=  max_y; n++)
			for (n_mode = 0; n_mode < Num_Prop; n_mode++)
				y_prop[n][n_mode] = 0.0;

		for (n_sph = 0; n_sph < num_obj; n_sph++)
		{
			objects[n_sph].p.x = 0.0;
			objects[n_sph].p.y = 0.0;
			objects[n_sph].p.z = 0.0;
			objects[n_sph].l.x = 0.0;
			objects[n_sph].l.y = 0.0;
			objects[n_sph].l.z = 0.0;
		}	
	}

	//--------------
	// Deinitialize 
	//--------------
	errcode = vslDeleteStream( &rngstream );
	CheckVslError( errcode );

	free(velcs_df[0][0]);
	free(velcs_df[0]);
	free(velcs_df);
	free(node_map[0]);
	free(node_map);
	free(objects);
	free(spheres);
	free(monomers);
	free(faces);
  free(forceDen[0]);  // Modification 20170420
  free(forceDen);
	fprintf (stdout, "End driver: proc #:%d, task %d\n", n_proc, task);
	fflush  (stdout);
}
