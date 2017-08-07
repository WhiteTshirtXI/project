/* Following Ahlrichs and Dunweg, modeling polymer segments as point forces in the fluid */
/* Coupled with LBE code to add hydrodynamic forces */


#include "header.h"

void spring_force(struct chain_param *chain_pm, struct monomer *mon);
void ev_force(Float radius, struct chain_param *chain_pm, struct monomer *mon);
void ev_force_nlist(Float radius, struct chain_param *chain_pm, struct monomer *monomers);
void wall_force(struct chain_param *chain_pm, struct monomer *mon);
void hi_force(struct monomer *mon, Float ***velcs_df, double fric, double tau, double dt, int num_beads,Float mon_mass, VSLStreamStatePtr);
void fluid_vel(int, double **, int x[], Float ***, struct monomer *, double *);
void fluc_force(Float *force, int n, VSLStreamStatePtr rngstream);

void get_forces(struct chain_param *chain_pm, struct monomer *monomers, Float ***velcs_df, int n_step, VSLStreamStatePtr rngstream)
{
  extern int n_proc, num_proc;
  extern int num_x;
  extern int max_x, max_y, max_z;
  extern int wall_flag;

  int mon_proc;
  int num_beads, Nbead;
  int i, j;
  int n1, n2;
  int s_type = 1;    /* spring type : 0 = FENE, 1=WLC */
  int ev_type = 1;   /* ev type: 0=HS, 1=WCA, 2=Gaussian */
  int maxsize[DIMS];
  int max_N_flag=0;
  double send_buffer[2*DIMS+1];
  char filename[80];
  FILE *stream;

  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  num_beads = chain_pm->num_beads;
  Nbead = chain_pm->num_beads / chain_pm->Nchain;

  for(j=0; j<DIMS; j++) 
    monomers[0].force[j] = monomers[0].f_ext[j];
  for(n1=1; n1<num_beads; n1++)
    for(j=0; j<DIMS; j++) 
      monomers[n1].force[j] = 0.0;

  /* Calculate spring and exc. vol. forces if # of beads > 1 */
  if(num_beads > 1) {
    spring_force(chain_pm, monomers);  /* 0 = FENE, 1=WLC */

    /*excluded volume forces */
    if(num_beads < 200)
      ev_force(monomers[0].radius, chain_pm, monomers);
    else {
      for(n1=0; n1<num_beads; n1++)
        if(monomers[n1].list[0] >= MAX_N-1) max_N_flag = 1;

      if(max_N_flag == 1)
        ev_force(monomers[0].radius, chain_pm, monomers);
      else ev_force_nlist(monomers[0].radius, chain_pm, monomers);
    }
  }

  for(n1=0; n1<num_beads; n1++)
    for(j=0; j<DIMS; j++) {
      monomers[n1].force[j] *= chain_pm->monmass;
      monomers[n1].f_int[j] = monomers[n1].force[j];
    }

  /*hydrodynamic and fluctuation forces */
  hi_force(monomers, velcs_df, chain_pm->fric, chain_pm->MD_steps, chain_pm->dt, num_beads, chain_pm->monmass, rngstream);
  
  if(wall_flag > 0)
    wall_force(chain_pm, monomers);

  /* send forces and fluid_vel back to proc 0 then broadcast from 0 */

  for(n1=0; n1< num_beads; n1++) { 
    mon_proc = (int)(monomers[n1].pos[0]/num_x);
    send_buffer[0]=n1;
    for(i=0; i<DIMS; i++) {
      send_buffer[i+1]=monomers[n1].force[i];
      send_buffer[i+DIMS+1]=monomers[n1].fluid_vel[i];
    }

    vector_copy(send_buffer, 2*DIMS+1, n_proc, mon_proc, 0);

    for(i=0; i<DIMS; i++) {
      monomers[(int)send_buffer[0]].force[i]=send_buffer[i+1];
      monomers[(int)send_buffer[0]].fluid_vel[i]=send_buffer[i+1+DIMS];
    }

  }
  sync_procs();

  for(n1=0; n1< num_beads; n1++) {
    broad_cast(monomers[n1].force, DIMS, 0);
    broad_cast(monomers[n1].fluid_vel, DIMS, 0);
  }
}


/* Calculate spring force : type = 0 FENE, 1 WLC */
void spring_force(struct chain_param *chain_pm, struct monomer *monomers)
{
  extern int wall_flag;
  extern int n_proc;
  extern int max_x, max_y, max_z;
  extern int num_x;
  long n1, n2;
  int i,j;
  int num_beads = chain_pm->num_beads;
  Float r1;
  Float slope, intercept;
  Float mag;
  Float Nks = chain_pm->nks*chain_pm->sigma_k;
  Float q12[DIMS];
  Float q12mag, q0;
  Float force[DIMS];
  int maxsize[DIMS];
  FILE *stream;

  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  for(n1=0; n1<num_beads-1; n1++) {
    n2=n1+1;
      
    if(monomers[n1].chain_id == monomers[n2].chain_id) {
      q12mag=0.0;

      for(j=0; j<DIMS; j++) 
	q12[j]=monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
	
      if(wall_flag < 3)
	q12[0]=n_image(q12[0], maxsize[0]);
      if(wall_flag < 2)
	q12[2]=n_image(q12[2], maxsize[2]);
      if(wall_flag < 1)
	q12[1]=n_image(q12[1], maxsize[1]);

      for(j=0; j<DIMS; j++)
	q12mag+=q12[j]*q12[j];
	
      q12mag = sqrt(q12mag);
	
      for(i=0; i<DIMS; i++)
	q12[i] /= q12mag;
	
      /* fene spring */
      if(chain_pm->spring == 0) {
	q0 = chain_pm->Q_fene*(monomers[n1].radius+monomers[n2].radius);
	  
	if(q12mag < q0) 
	  mag = chain_pm->kT*chain_pm->H_fene*q12mag/(1.0 - (q12mag/q0)*(q12mag/q0));
	else {
	  stream = fopen("output.err", "w");
	  fprintf(stream, "Spring is broken, %d %le %le\n", n1, q12mag, q0);
	  fprintf(stream, "mon %d pos=(%le %le %le), mon %d pos = (%le %le %le)\n", n1, monomers[n1].pos_pbc[0], monomers[n1].pos_pbc[1], monomers[n1].pos_pbc[2], n2, monomers[n2].pos_pbc[0], monomers[n2].pos_pbc[1], monomers[n2].pos_pbc[2]);
	  fprintf(stream, "mon %d force=(%le %le %le), mon %d force = (%le %le %le)\n", n1, monomers[n1].force[0], monomers[n1].force[1], monomers[n1].force[2], n2, monomers[n2].force[0], monomers[n2].force[1], monomers[n2].force[2]);
	  fclose(stream);
	  exit(11);
	}
      }
      
      /* worm like spring */
      else if (chain_pm->spring == 1) {
	r1 = q12mag/Nks;
	if(r1 < 1.0)
	  mag = 2.0*chain_pm->kT/chain_pm->sigma_k*(0.25/((1.0-r1)*(1.0-r1))-0.25+r1);
	else {	  
	  if(n_proc == 0) {
	    stream = fopen("output.err", "w");
	    fprintf(stream, "Spring is broken, %d %le\n", n1, q12mag);
	    fprintf(stream, "mon %d pos=(%le %le %le), mon %d pos = (%le %le %le)\n", n1, monomers[n1].pos_pbc[0], monomers[n1].pos_pbc[1], monomers[n1].pos_pbc[2], n2, monomers[n2].pos_pbc[0], monomers[n2].pos_pbc[1], monomers[n2].pos_pbc[2]);
	    fprintf(stream, "mon %d force=(%le %le %le), mon %d force = (%le %le %le)\n", n1, monomers[n1].force[0], monomers[n1].force[1], monomers[n1].force[2], n2, monomers[n2].force[0], monomers[n2].force[1], monomers[n2].force[2]);
	    fclose(stream);
	    exit(11);
	  }
	}
      }
	
      /* harmonic spring */
      else if (chain_pm->spring == 2)
	mag = 2.0*chain_pm->kT*chain_pm->H_fene*(q12mag-1.0);
	
      for(i=0; i<DIMS; i++) {
	force[i]=q12[i]*mag;
	  
	monomers[n1].force[i]-=force[i];
	monomers[n2].force[i]+=force[i];
      }
    }
  }
}

/* Exc. Vol forces type = 0 HS, 1 LJ, 2 Gaussian */
void ev_force(Float radius, struct chain_param *chain_pm, struct monomer *monomers)
{
  extern int wall_flag;
  extern int n_proc;
  extern int max_x, max_y, max_z;
  extern int num_x;
  extern double tau;
  int i, j;
  long n1, n2;
  int num_beads = chain_pm->num_beads;
  Float sigma = 2.0*radius;
  Float cutoff = chain_pm->evcutoff*radius;
  Float aterm;
  Float q12[DIMS];
  Float q12mag;
  Float eps, Ecut;
  Float force[DIMS];
  int maxsize[DIMS];
  FILE *stream;
  double eta = (tau - 0.5)/3.0;
  double Rg = chain_pm->Ss;
  double Nks =chain_pm->nks;
  double v = chain_pm->sigma_k*chain_pm->sigma_k*chain_pm->sigma_k;
  int Nbead = chain_pm->num_beads/chain_pm->Nchain;
  double aev_alpha = 4.0/3.0*Rg*Rg;
  double aev = Nks*Nks*(1.0/((Pi*aev_alpha)*sqrt(Pi*aev_alpha)));
		
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  cutoff = 1.12*sigma;
  eps = chain_pm->kT/sigma;

  for(i=0; i<DIMS; i++)
    force[i]=0.0;

  for(n1=0; n1 < num_beads-1; n1++) {
    for(n2=n1+1; n2 < num_beads; n2++) {

      q12mag=0.0;
      for(j=0; j<DIMS; j++) 
	q12[j]=monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
	
      if(wall_flag < 3)
	q12[0]=n_image(q12[0], maxsize[0]);
      if(wall_flag < 2)
	q12[2]=n_image(q12[2], maxsize[2]);
      if(wall_flag < 1)
	q12[1]=n_image(q12[1], maxsize[1]);

      for(j=0; j<DIMS; j++) 
	q12mag+=q12[j]*q12[j];
	
      q12mag = sqrt(q12mag);
	
      if(chain_pm->ev_type ==0) {
	for(i=0; i<DIMS; i++)
	  force[j]=0.0;
      }
	
      /* WCA potential */
      else if(chain_pm->ev_type ==1) {
	if(q12mag > cutoff) 
	  for(i=0; i<DIMS; i++)
	    force[i]=0.0;
	else if(q12mag < cutoff) {
	  aterm = (sigma/q12mag);
	  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
	  
	  for(i=0; i<DIMS; i++) 
	    force[i]=-24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
	}
      }
      
      else if(chain_pm->ev_type ==2) {
	if(q12mag > chain_pm->evcutoff*Rg)
	  for(i=0; i<DIMS; i++)
	    force[i]=0.0;
	else {
	  aterm = aev/aev_alpha * exp(-q12mag*q12mag/aev_alpha);
	    
	  for(i=0; i<DIMS; i++) 
	    force[i]= (-v*chain_pm->kT*aterm*q12[i]);
	}	
      }
	
      for(i=0; i<DIMS; i++) {
	monomers[n1].force[i] -= force[i];
	monomers[n2].force[i] += force[i];
      }
    }
  }
}


/* Exc. Vol forces type = 0 HS, 1 LJ, 2 Gaussian */
void ev_force_nlist(Float radius, struct chain_param *chain_pm, struct monomer *monomers)
{
  extern int wall_flag;
  extern int n_proc;
  extern int max_x, max_y, max_z;
  extern int num_x;
  extern double tau;
  int i, j;
  long n1, n2;
  int num_beads = chain_pm->num_beads;
  Float sigma = 2.0*radius;
  Float cutoff = chain_pm->evcutoff*radius;
  Float aterm;
  Float q12[DIMS];
  Float q12mag;
  Float eps, Ecut;
  Float force[DIMS];
  int maxsize[DIMS];
  FILE *stream;
  double eta = (tau - 0.5)/3.0;
  double Rg = chain_pm->Ss;
  double Nks =chain_pm->nks;
  double v = chain_pm->sigma_k*chain_pm->sigma_k*chain_pm->sigma_k;
  int Nbead = chain_pm->num_beads/chain_pm->Nchain;
  double aev_alpha = 4.0/3.0*Rg*Rg;
  double aev = Nks*Nks*(1.0/((Pi*aev_alpha)*sqrt(Pi*aev_alpha)));
		
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  cutoff = 1.12*sigma;
  eps = chain_pm->kT/sigma;
  aterm = (sigma/cutoff);
  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
  Ecut=eps*(2.0*(aterm*aterm)-aterm);

  for(i=0; i<DIMS; i++)
    force[i]=0.0;

  for(n1=0; n1 < num_beads; n1++) {
    for(n2=1; n2 <= monomers[n1].list[0]; n2++) {
      q12mag=0.0;

      for(j=0; j<DIMS; j++) 
	q12[j]=monomers[n1].pos_pbc[j] - monomers[monomers[n1].list[n2]].pos_pbc[j];
	
      if(wall_flag < 3)
	q12[0]=n_image(q12[0], maxsize[0]);
      if(wall_flag < 2)
	q12[2]=n_image(q12[2], maxsize[2]);
      if(wall_flag < 1)
	q12[1]=n_image(q12[1], maxsize[1]);

      for(j=0; j<DIMS; j++) 
	q12mag+=q12[j]*q12[j];
	
      q12mag = sqrt(q12mag);

//      if(chain_pm->spring == 2)
//	if(monomers[n1].list[n2]==n1+1 && monomers[n1].chain_id == monomers[monomers[n1].list[n2]].chain_id)
//	  continue;
	
      if(chain_pm->ev_type ==0) {
	for(i=0; i<DIMS; i++)
	  force[j]=0.0;
      }
	
      /* WCA potential */
      else if(chain_pm->ev_type ==1) {
	if(q12mag > cutoff) 
	  for(i=0; i<DIMS; i++)
	    force[i]=0.0;
	else if(q12mag < cutoff) {
	  aterm = (sigma/q12mag);
	  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
	      
	  for(i=0; i<DIMS; i++) 
	    force[i]=-24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
	}
      }
	
      else if(chain_pm->ev_type ==2) {
	aterm = aev/aev_alpha * exp(-q12mag*q12mag/aev_alpha);
	for(i=0; i<DIMS; i++) 
	  force[i]= (-v*chain_pm->kT*aterm*q12[i]);
      }
	
      
      for(i=0; i<DIMS; i++) {
	/*
	  if(fabs(force[i]) > 100)
	  printf("excluded volume force too large, (%d %d), (%le %le %le), (%le %le %le)\n", n1, n2, monomers[n1].pos[0], monomers[n1].pos[1], monomers[n1].pos[2], monomers[monomers[n1].list[n2]].pos[0], monomers[monomers[n1].list[n2]].pos[1], monomers[monomers[n1].list[n2]].pos[2]);
	*/

	monomers[n1].force[i] -= force[i];
      /*monomers[monomers[n1].list[n2]].force[i] += force[i]; */
      }
    }
  }
}

/* Calculate spring force : type = 0 FENE, 1 WLC */
void wall_force(struct chain_param *chain_pm, struct monomer *monomers)
{
  extern int n_proc;
  extern int num_x;
  extern int max_x, max_y, max_z;
  extern int wall_flag;
  long n1, n2;
  int i,j;
  int num_beads = chain_pm->num_beads;
  Float r1;
  Float mag;
  Float wally, wallz;
  Float wall_aev;
  Float force[DIMS];
  double maxsize[DIMS], monpos[DIMS];
  FILE *stream;

  maxsize[0]=(double)max_x;
  maxsize[1]=(double)max_y;
  maxsize[2]=(double)max_z;
  
  if(monomers[0].radius*3 > 1.0) {
    //    wally = (maxsize[1]-1.0)/2.0-monomers[0].radius;
    //    wallz = (maxsize[2]-1.0)/2.0-monomers[0].radius;
    wally = (maxsize[1]-1.0)/2.0-1.2;
    wallz = (maxsize[2]-1.0)/2.0-1.2;
    wall_aev = 25.0*chain_pm->kT/(chain_pm->sigma_k*chain_pm->sigma_k*chain_pm->sigma_k);
  }
  else {
    wally = (maxsize[1]-1.0)/2.0-1.0;
    wallz = (maxsize[2]-1.0)/2.0-1.0;
    wall_aev = 25.0*chain_pm->kT/(chain_pm->sigma_k*chain_pm->sigma_k*chain_pm->sigma_k);
  }

  for(n1=0; n1<num_beads; n1++) {
    if((int)(monomers[n1].pos[0]/num_x) == n_proc)
      monpos[0] = monomers[n1].pos[0] - ((maxsize[0]-1.0)/2.0+1.0);
    for(j=1; j<DIMS; j++)
      monpos[j] = monomers[n1].pos[j] - ((maxsize[j]-1.0)/2.0);

    if(wall_flag >= 1) {
      if(monpos[1] < -wally)
	force[1] = wall_aev * (monpos[1]+wally)*(monpos[1]+wally);
      else if(monpos[1] > wally)
	force[1] = -wall_aev * (monpos[1]-wally)*(monpos[1]-wally);
      else
	force[1]=0.0;
      
      monomers[n1].force[1] += force[1];

      if(fabs(monpos[1]) > (maxsize[1]-1.0)/2.0) {
	printf(" bead %d out of ybounds, pos=%le %le %le \n", n1, monomers[n1].pos[0], monomers[n1].pos[1], monomers[n1].pos[2]);
	printf(" fluid velocity = %le %le %le  vel=%le %le %le\n", monomers[n1].fluid_vel[0], monomers[n1].fluid_vel[1], monomers[n1].fluid_vel[2], monomers[n1].vel[0], monomers[n1].vel[1], monomers[n1].vel[2]);
	printf(" force = %le %le %le internal force=%le %le %le\n", monomers[n1].force[0], monomers[n1].force[1], monomers[n1].force[2], monomers[n1].f_int[0], monomers[n1].f_int[1], monomers[n1].f_int[2]);
	printf(" wall force = %le \n", force[1]);
	exit(18);
      }
    }

    if(wall_flag == 2) {
      if(monpos[1] < -wally)
	force[1] = wall_aev * (monpos[1]+wally)*(monpos[1]+wally);
      else if(monpos[1] > wally)
	force[1] = -wall_aev * (monpos[1]-wally)*(monpos[1]-wally);
      else
	force[1]=0.0;

      if(monpos[2] > wallz)
	force[2] = -wall_aev * (monpos[2]-wallz)*(monpos[2]-wallz);
      else if(monpos[2] < -wallz)
	force[2] = wall_aev * (monpos[2]+wallz)*(monpos[2]+wallz);
      else
	force[2] = 0.0;
      
      monomers[n1].force[1] += force[1];
      monomers[n1].force[2] += force[2];

      if(fabs(monpos[1]) > (maxsize[1]-1.0)/2.0 || fabs(monpos[2]) > (maxsize[2]-1.0)/2.0) {
	printf(" bead %d out of yzbounds, pos=%le %le %le \n", n1, monomers[n1].pos[0], monomers[n1].pos[1], monomers[n1].pos[2]);
	printf(" fluid velocity = %le %le %le  vel=%le %le %le\n", monomers[n1].fluid_vel[0], monomers[n1].fluid_vel[1], monomers[n1].fluid_vel[2], monomers[n1].vel[0], monomers[n1].vel[1], monomers[n1].vel[2]);
	printf(" force = %le %le %le internal force=%le %le %le\n", monomers[n1].force[0], monomers[n1].force[1], monomers[n1].force[2], monomers[n1].f_int[0], monomers[n1].f_int[1], monomers[n1].f_int[2]);
	printf(" wall force = %le %le\n", force[1], force[2]);
	exit(18);
      }
    }
  }

}


/* Hydrodynamic forces */
void  hi_force(struct monomer *mon, Float ***velcs_df, double fric, double tau, double dt, int num_beads, Float mon_mass, VSLStreamStatePtr rngstream)
{
  extern struct vector f_ext;
  extern int n_proc, num_proc, num_x;
  extern int max_x, max_y, max_z;
  extern int wall_flag;
  extern int add_noise;
  extern int backflow_flag;
  static int MD_tau=0;

  int mon_proc;
  int i, j, k, q, d;
  long n, nxy, n1;
  int x[DIMS], x1[27][DIMS];
  int qxp, qyp, qzp;
  int qxm, qym, qzm;
  int maxsize[DIMS];
  double x2[DIMS], x3[DIMS], x4[DIMS];
  double **vel, *tmp_pp;
  double weight[27], vel_node[27][DIMS], rho_node[27];
  double temp[DIMS], temp1[DIMS];
  double radius = mon[0].radius;
  double sigma = 2.0*radius;
  double vol;
  double dmass, dmom[3];
  Float dmomentum[Num_Dir], dj_temp[DIMS];
  Float *force, *fforce;

  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  vol=max_x*max_y*max_z;
 
  force=(Float *)calloc(DIMS, sizeof(Float));
  if(force == 0) fatal_err("cannot allocate force", -1);

  fforce=(Float *)calloc(DIMS*num_beads, sizeof(Float));
  if(fforce == 0) fatal_err("cannot allocate force", -1);

  vel = (double **)calloc (num_beads, sizeof(*vel));
  if (vel == 0)  fatal_err ("cannot allocate vel[][]", -1);
  
  tmp_pp = (double *)calloc (num_beads*DIMS, sizeof(*tmp_pp));
  if (tmp_pp == 0)  fatal_err ("cannot allocate vel[][]", -1);
  for (n1 = 0; n1 < num_beads; n1++)                          /* Assign pointers to pointers */
    {
      vel[n1] = tmp_pp;
      tmp_pp += DIMS;
    }

  /* calculate drag force */
  /* calculate fluid velocity at the monomer coords. */
  for(n1=0; n1< num_beads; n1++) {
    for(d=0; d<DIMS; d++) {
      x[d] = mon[n1].pos[d];
      if(mon[n1].pos[d]-x[d] > 0.5)
	x[d]++;
    }

    fluid_vel(n1, vel, x, velcs_df, mon, weight);
    
    //    if(MD_tau == tau)
    for(d=0; d<DIMS; d++) 
      mon[n1].fluid_vel[d]=vel[n1][d];

    for(d=0; d<DIMS; d++) 
      mon[n1].fricforce[d]=-fric*(mon[n1].vel[d]-vel[n1][d]);
      //mon[n1].fricforce[d]=-fric*(mon[n1].vel[d]);
  }

  /* Generate the random force */
  if(add_noise >= 1 && add_noise !=3)
    fluc_force(fforce, num_beads*DIMS, rngstream);

  /* redistribute momentum to the surrounding lattice */
  n1=0;
  while(n1<num_beads) {
    for(d=0; d<DIMS; d++) {
      x[d] = (int)mon[n1].pos[d];
      if(mon[n1].pos[d]-x[d] > 0.5)
	x[d]=x[d]+1;
    }

    mon_proc = box(x[0],maxsize[0])/num_x;
    if(n_proc == mon_proc) {
      for(i=0; i<DIMS; i++) {
	force[i] = Rho_Fl*mon[n1].fricforce[i];
	force[i] += fforce[n1*DIMS+i];
      }
      
      for(i=0; i<DIMS; i++) {
	mon[n1].f_fluc[i]=fforce[n1*DIMS+i];
	mon[n1].force[i] +=force[i];
      }
      
      dmass = 0.0;
      for(i=0; i<DIMS; i++) {
	dmom[i] = 0.0;
	//	dj_temp[i]=-force[i]/(tau*mon_mass*CS2*Rho_Fl);
	dj_temp[i]=-force[i]*dt/(mon_mass*CS2*Rho_Fl);
	//dj_temp[i]=-force[i]*dt/(CS2*Rho_Fl);
      }
      
      // printf("vel=%le fluidvel=%le force=%le friction=%le\n", mon[n1].vel[0], mon[n1].fluid_vel[0], mon[n1].force[0], force[0]);
      // printf("%le %le %le %le\n", mon[n1].vel[0], mon[n1].fluid_vel[0], mon[n1].force[0], force[0]);
      
      for(q=0; q<Num_Dir; q++) 
	dmomentum[q] = fac[q]*(c_x[q] * dj_temp[0] + c_y[q] * dj_temp[1] + c_z[q] * dj_temp[2]);
      
      /* Use trilinear interpolation to get the fluid velocity at the monomer position */
      for(i=0; i<=2; i++) { 
	for(j=0; j<=2; j++) {
	  for(k=0; k<=2; k++) {
	    n = 9*i+3*j+k;
	    x1[n][2] = box(x[2]+(k-1), maxsize[2])+1;
	    x1[n][1] = box(x[1]+(j-1), maxsize[1])+1;
	    x1[n][0] = box((x[0]+(i-1)), maxsize[0])%num_x+1.0;
	    
	    nxy = x1[n][0]*(maxsize[1]+2)+(x1[n][1]);
	    
	    for(q=0; q<Num_Dir; q++) 
	      velcs_df[nxy][q][x1[n][2]] += dmomentum[q]*weight[n]; 
       
/* 	    for(d=0; d<DIMS; d++) */
/* 	      vel_node[n1][d] = 0.0; */
/* 	    for(q=0; q<Num_Dir; q++) { */
/* 	      vel_node[n1][0] += c_x[q]*velcs_df[nxy][q][x1[2]]; */
/* 	      vel_node[n1][1] += c_y[q]*velcs_df[nxy][q][x1[2]]; */
/* 	      vel_node[n1][2] += c_z[q]*velcs_df[nxy][q][x1[2]]; */

/*       	      dmass += velcs_df[nxy][q][x1[2]]; */
/* 	    } */

/* 	    dmom[0] += vel_node[n1][0]; */
/* 	    dmom[1] += vel_node[n1][1]; */
/* 	    dmom[2] += vel_node[n1][2]; */
	  
/* 	    printf("node (%d %d %d) vel (%le %le %le)\n", x1[0], x1[1], x1[2], vel_node[n1][0], vel_node[n1][1], vel_node[n1][2]); */
  	  } 
 	} 
      }
      /*       printf("dmass=%le dmom=(%le %le %le), djtemp=(%le %le %le)\n", dmass, dmom[0], dmom[1], dmom[2], dj_temp[0], dj_temp[1], dj_temp[2]); */
    }
    n1++;
  }

  if(MD_tau == tau)
    MD_tau =1;
  else 
    MD_tau++;

  free(force);
  free(fforce);
  free(vel[0]);
  free(vel);
}

void fluid_vel(int n1, double **vel, int x[DIMS], Float ***velcs_df, struct monomer *mon, double weight[])
{
  extern int wall_flag;
  extern int max_x, max_y, max_z;
  extern int num_x;

  int i,j,k;
  int d,q,n,nxy;
  int x1[27][DIMS], maxsize[DIMS];
  double x2[DIMS], x3[DIMS], x4[DIMS];
  double vel_node[27][DIMS], rho_node[27];
  double temp[DIMS];
  double rho;
  
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  for(i=0; i<27; i++)
    rho_node[i]=0.0;

  for(i=0; i<=2; i++) {
    for(j=0; j<=2; j++) {
      for(k=0; k<=2; k++) {
  	n = 9*i+3*j+k;
    	x1[n][2] = box(x[2]+(k-1), maxsize[2])+1;
	x1[n][1] = box(x[1]+(j-1), maxsize[1])+1;
	x1[n][0] = box((x[0]+(i-1)), maxsize[0])%num_x+1.0;
      
	nxy = x1[n][0]*(maxsize[1]+2)+(x1[n][1]);
	
	rho_node[n]=0.0;
	for(d=0; d<DIMS; d++)
	  vel_node[n][d]=0.0;
	
	for(q=0; q<Num_Dir; q++) { 
	  rho_node[n]   +=velcs_df[nxy][q][x1[n][2]];
	  vel_node[n][0]+=c_x[q]*velcs_df[nxy][q][x1[n][2]];
	  vel_node[n][1]+=c_y[q]*velcs_df[nxy][q][x1[n][2]];
	  vel_node[n][2]+=c_z[q]*velcs_df[nxy][q][x1[n][2]];
	}
	
	// 	printf("pre node rho %le (%d %d %d) vel (%le %le %le)\n", rho_node[n], x1[n][0], x1[n][1], x1[n][2], vel_node[n][0], vel_node[n][1], vel_node[n][2]); 
      }
    }
  }

  /* Use trilinear interpolation to get the fluid velocity at the monomer position */
  for(d=0; d<DIMS; d++) {
    vel[n1][d] = 0.0;
    x2[d]=mon[n1].pos[d]-(double)x[d];   /* translate the coordinates of the monomer */
    x3[d]=1.0+x2[d];
    x4[d]=1.0-x2[d];

    x2[d]=fabs(x2[d]);
  }

  for(i=0; i<=2; i++) { 
    if(i==0)
      temp[0]=(5.0-3*x3[0]-sqrt(6*x3[0]-3*x3[0]*x3[0]-2.0))/6.0;
    else if(i==1) 
      temp[0]=(1.0+sqrt(1-3*x2[0]*x2[0]))/3.0;
    else if(i==2)
      temp[0]=(5.0-3*x4[0]-sqrt(6*x4[0]-3*x4[0]*x4[0]-2.0))/6.0;
    
    for(j=0; j<=2; j++) {
      if(j==0)
	temp[1]=(5.0-3*x3[1]-sqrt(6*x3[1]-3*x3[1]*x3[1]-2.0))/6.0;
      else if(j==1) 
	temp[1]=(1.0+sqrt(1-3*x2[1]*x2[1]))/3.0;
      else if(j==2)
	temp[1]=(5.0-3*x4[1]-sqrt(6*x4[1]-3*x4[1]*x4[1]-2.0))/6.0;
	
      for(k=0; k<=2; k++) {
	if(k==0)
	  temp[2]=(5.0-3*x3[2]-sqrt(6*x3[2]-3*x3[2]*x3[2]-2.0))/6.0;
	else if(k==1) 
	  temp[2]=(1.0+sqrt(1-3*x2[2]*x2[2]))/3.0;
	else if(k==2)
	  temp[2]=(5.0-3*x4[2]-sqrt(6*x4[2]-3*x4[2]*x4[2]-2.0))/6.0;
	
	n=9*i+3*j+k;
	weight[n] = temp[0]*temp[1]*temp[2];
	for(d=0; d<DIMS; d++) {
	  vel[n1][d]+=vel_node[n][d]*weight[n];
	}
      }
    }
  }
  
  for(d=0; d<DIMS; d++)    
    vel[n1][d] /= Rho_Fl;
}

void fluc_force(Float *force, int n, VSLStreamStatePtr rngstream)
{
  extern Float mon_stdev;

  int errcode;

  /* Calculate the fluctuations */  
  errcode=vdRngGaussian( METHOD, rngstream, n, force, 0.0, mon_stdev);
  CheckVslError(errcode);

  //  printf("%le %le %le\n", f_fluc[0], f_fluc[1], f_fluc[2]);  
}
