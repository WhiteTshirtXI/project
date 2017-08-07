#include "header.h"

/* void fluc_vel(Float *force); */
void temp_rescale(struct monomer *mon, struct sphere_param *sphere_pm);

void verlet_update(struct monomer *mon, struct face *faces, struct sphere_param *sphere_pm, 
                   Float ***velcs_df, int n_step, VSLStreamStatePtr rngstream)
{
	extern int max_x, max_y, max_z;
	extern int wall_flag;
	double maxsize[3];
	long i,j;
	int d;
	int num_beads = sphere_pm->num_beads;
	int overlap;
	double dt = sphere_pm->dt;
	double h_dt = dt/2.0;
	double h_dt2= h_dt*h_dt;
	double trialpos[DIMS], trialvel[DIMS];
	double *fforce;
	double dr[DIMS];
	double wallx, wally, wallz;
	double mon_mass = Rho_Fl*sphere_pm->monmass;
	double x = sphere_pm->fric*dt/sphere_pm->monmass;
	double midpt[DIMS];
	FILE *stream;
  extern char *work_dir;  // Modification 20170321
  char fileName[100];
	maxsize[0]=(double)max_x;
	maxsize[1]=(double)max_y;
	maxsize[2]=(double)max_z;

	fforce=(Float *)calloc(DIMS, sizeof(Float));
	if(fforce == 0) fatal_err("cannot allocate force", -1);

	for(j=0; j<DIMS; j++)
		midpt[j]=(maxsize[j]-1.0)/2.0;

	wallx = (maxsize[0]+1.0)/2.0-1.5*mon[0].radius;
	wally = (maxsize[1])/2.0-1.5*mon[0].radius;
	wallz = (maxsize[2])/2.0-1.5*mon[0].radius;

	if(sphere_pm->verlet == 0) /* no position update */
	{
		get_forces(sphere_pm, mon, faces, velcs_df, n_step, rngstream);

		for(i=0;i<num_beads;i++) {
			for(d=0; d<DIMS; d++) {
				mon[i].vel[d] += mon[i].force[d]/mon_mass*dt;
				dr[d] = dt*trialvel[d];
				mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
			}
		}
	} 
	else if(sphere_pm->verlet==1)  /* velocity verlet update */
	{
		for(i=0;i < num_beads;i++)	
		{
			for(d=0; d < DIMS; d++) {
				dr[d] = dt * mon[i].vel[d] + 2.0 * h_dt2 * mon[i].force0[d]/mon_mass;
				mon[i].pos_pbc[d] = mon[i].pos_pbc[d] + dr[d];
				mon[i].pos[d] = box(mon[i].pos_pbc[d], maxsize[d]);

        //if(isnan(mon[i].pos_pbc[d]))
        //{
        //  sprintf(fileName,"%s/data/nan_warning.dat",work_dir);
        //  stream=fopen(fileName,"w");
        //  fprintf(stream,"step = %d  momomer %d\n",n_step,i);
        //  fprintf(stream,"spring force = (%f, %f, %f)\n",mon[i].force_spring[0], 
        //          mon[i].force_spring[1],mon[i].force_spring[2]);
        //  fprintf(stream,"bending force = (%f, %f, %f)\n",mon[i].force_bend[0],
        //          mon[i].force_bend[1],mon[i].force_bend[2]);
        //  fprintf(stream,"wall force = (%f, %f, %f)\n",mon[i].force_wall[0],
        //          mon[i].force_wall[1],mon[i].force_wall[2]);
        //  fprintf(stream,"drag foce = (%f, %f, %f)\n",mon[i].drag[0],
        //          mon[i].drag[1],mon[i].drag[2]);
        //  fprintf(stream,"inter-particle foce = (%f, %f, %f)\n",mon[i].force_pp[0],
        //          mon[i].force_pp[1],mon[i].force_pp[2]);
        //  fprintf(stream,"intra-particle nonbonded foce = (%f, %f, %f)\n",
        //          mon[i].force_nonbonded_intra[0], mon[i].force_nonbonded_intra[1],
        //          mon[i].force_nonbonded_intra[2]);
        //  fclose(stream);
        //  exit(18);
        //}
			}
		}
   
    //if(n_step > 11999) 
    //  ArtificialShift(sphere_pm, mon);  // Modification 20170323
    // Modification 20170423
		get_forces_hi(sphere_pm, mon, faces, velcs_df, n_step, rngstream);

		for(i=0;i < num_beads;i++) 
		{
			for(d=0; d < DIMS; d++) {
				mon[i].vel[d] += h_dt * (mon[i].force0[d]+mon[i].force[d])/mon_mass;
				mon[i].force0[d]=mon[i].force[d];
			}
		}
	}
	else if(sphere_pm->verlet == 2) /* explicit 1st order */
	{
		get_forces(sphere_pm, mon, faces, velcs_df, n_step, rngstream);

		for(i=0;i<num_beads;i++) {
			for(d=0; d<DIMS; d++) {
				trialvel[d] = (mon[i].vel[d]+dt*mon[i].force[d]/mon_mass);
				dr[d] = dt*(trialvel[d]+mon[i].vel[d])/2.0;
				trialpos[d]=mon[i].pos_pbc[d]+dr[d];
			}

			for(d=0; d<DIMS; d++) {
				mon[i].vel[d]=trialvel[d];
				mon[i].pos_pbc[d]=trialpos[d];
				mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
			}
		}
	}
	else if(sphere_pm->verlet == 3) /* implicit 1st order */
	{
		for(i=0;i<num_beads;i++) {
			for(d=0; d<DIMS; d++) {
				mon[i].pos_tmp[d]=mon[i].pos_pbc[d];
				mon[i].pos_pbc[d] +=mon[i].vel[d]*dt;
				mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
			}
		}

		get_forces(sphere_pm, mon, faces, velcs_df, n_step, rngstream);

		for(i=0;i<num_beads;i++) {
			for(d=0; d<DIMS; d++) {
				trialvel[d] = (mon[i].vel[d]+dt*(mon[i].force[d])/mon_mass);
				dr[d] = dt*(trialvel[d]+mon[i].vel[d])/2.0;
				mon[i].vel[d]=trialvel[d];
				mon[i].pos_pbc[d]=mon[i].pos_tmp[d]+dr[d];
				mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
			}

			/*
				 for(d=0; d<DIMS; d++) {
				 trialvel[d] = (mon[i].vel[d]+dt*(mon[i].f_int[d]+mon[i].f_fluc[d])/mon_mass)/(1+x);
				 dr[d] = dt*(trialvel[d]+mon[i].vel[d])/2.0;
				 mon[i].vel[d]=trialvel[d];
				 mon[i].pos_pbc[d]=mon[i].pos_tmp[d]+dr[d];
				 mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
				 }

				 for(d=0; d<DIMS; d++) {
				 mon[i].vel[d]=(1+x)/(1+2*x)*(mon[i].vel[d]+mon[i].fluid_vel[d]*x/(1+x)+
				 dt*(mon[i].f_int[d]/mon_mass+mon[i].f_fluc[d]/mon_mass/(1+x)));
				 dr[d] = dt*mon[i].vel[d];
				 mon[i].pos_pbc[d]=mon[i].pos_tmp[d]+dr[d];
				 mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
				 }
			 */
		}
	}

	free(fforce);
}

void temp_rescale(struct monomer *mon, struct sphere_param *sphere_pm)
{
	int i, d;
	double kT_ob; 
	double scale;
	double avg_vel2[DIMS];

	for(d=0; d<DIMS; d++)
		avg_vel2[d] = 0.0;

	for(i=0; i<sphere_pm->num_beads; i++) {
		for(d=0; d<DIMS; d++)
			avg_vel2[d] += mon[i].vel[d]*mon[i].vel[d];
	}

	kT_ob = 0.0;
	for(d=0; d<DIMS; d++)
		kT_ob += avg_vel2[d];

	kT_ob /= sphere_pm->num_beads;

	if(kT_ob < 1e-10)
		scale = 0.0;
	else
		scale = sphere_pm->kT / (kT_ob*Rho_Fl/DIMS);

	sphere_pm->tempscale = scale;
}

/* void fluc_vel(Float *force) */
/* { */
/*   extern mt_state twister; */
/*   extern Float mon_stdev; */

/*   int d; */
/*   double f_fluc[DIMS]; */

/*   /\* Calculate the fluctuations *\/ */
/*   for(d=0; d<DIMS; d++) { */
/*     f_fluc[d] = rds_normal(&twister, 0.0e0, mon_stdev); */
/*     force[d] =f_fluc[d]/Rho_Fl; */
/*   } */

/*   //  printf("%le %le %le\n", f_fluc[0], f_fluc[1], f_fluc[2]);   */
/* } */

double adams_bashforth (double y, double derivOld, double derivNew, double dt) 
{
  double yNew;
  yNew = y + (1.5*derivNew - 0.5*derivOld) * dt;
  return yNew;
}

double euler_method (double y, double deriv, double dt) 
{
  double yNew;
  yNew = y + deriv * dt;
  return yNew;
}

void update_position (int numBead, double dt, struct monomer *mon) 
{
  extern int max_x, max_y, max_z;
  double maxsize[DIMS];
  maxsize[0] = max_x;
  maxsize[1] = max_y;
  maxsize[2] = max_z;
  for(int n=0; n < numBead; n++) 
  {
    if(mon[n].updatedFlag==FALSE)
    {
      mon[n].pos_pbc[0] = euler_method(mon[n].pos_pbc[0], mon[n].vel[0], dt);
      //mon[n].pos_pbc[d] = AdamsBashforth(mon[n].pos_pbc[d], mon[n].velOld[d], mon[n].vel[d], dt);
      mon[n].pos_pbc[1] = euler_method(mon[n].pos_pbc[1], mon[n].vel[1], dt);
      mon[n].pos_pbc[2] = euler_method(mon[n].pos_pbc[2], mon[n].vel[2], dt);
      mon[n].pos[0] = box(mon[n].pos_pbc[0], maxsize[0]);
      mon[n].pos[1] = box(mon[n].pos_pbc[1], maxsize[1]);
      mon[n].pos[2] = box(mon[n].pos_pbc[2], maxsize[2]);
      mon[n].vel_old[0] = mon[n].vel[0];
      mon[n].vel_old[1] = mon[n].vel[1];
      mon[n].vel_old[2] = mon[n].vel[2];
    }
    //mon[n].updatedFlag = FALSE;  // Modification 20170723 move to get_force.c
  }
}

void euler_update(int numBead, struct monomer *mon)
{
  // Modification 20170712
  extern double fictionalMass;
  extern int max_x, max_y, max_z;
  double maxsize[DIMS];
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  for(int n=0; n < numBead; n++) {
    mon[n].vel[0] += mon[n].force[0]/fictionalMass;
    mon[n].pos_pbc[0] += mon[n].vel[0];
    mon[n].pos[0] = box(mon[n].pos_pbc[0], maxsize[0]);
    mon[n].vel[1] += mon[n].force[1]/fictionalMass;
    mon[n].pos_pbc[1] += mon[n].vel[1];
    mon[n].pos[1] = box(mon[n].pos_pbc[1], maxsize[1]);
    mon[n].vel[2] += mon[n].force[2]/fictionalMass;
    mon[n].pos_pbc[2] += mon[n].vel[2];
    mon[n].pos[2] = box(mon[n].pos_pbc[2], maxsize[2]);
  }
}

