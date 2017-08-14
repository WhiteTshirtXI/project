/* Following Ahlrichs and Dunweg, modeling polymer segments as point forces in the fluid */
/* Coupled with LBE code to add hydrodynamic forces */

#include "header.h"
void spring_force(struct sphere_param *sphere_pm, struct monomer *mon, struct face *faces);
void ev_force(Float radius, struct sphere_param *sphere_pm, struct monomer *mon);
void ev_force_nlist(Float radius, struct sphere_param *sphere_pm, struct monomer *monomers);
void wall_force(struct sphere_param *sphere_pm, struct monomer *mon);
void hi_force(struct monomer *mon, Float ***velcs_df, double fric, double tau, double dt, 
              int num_beads,Float mon_mass, VSLStreamStatePtr);
void fluid_vel(int, double **, int x[], Float ***, struct monomer *);
void fluc_force(Float *force, int n, VSLStreamStatePtr rngstream);
void product(double a[DIMS], double b[DIMS], double c[DIMS]);
void ParticleStress(struct sphere_param *sphere_pm, struct monomer *monomers);
void zero_variable(int numBead,struct monomer *monomers);

void zero_variable(int numBead,struct monomer *monomers)
{
	for(int n1=0; n1 < numBead; n1++) 
  {
    monomers[n1].updatedFlag = FALSE;
		for(int j=0; j<DIMS; j++) { 
			monomers[n1].force[j]=0;
      monomers[n1].force_spring[j]=0.;
      monomers[n1].force_bend[j]=0.;
      monomers[n1].force_vol[j]=0.;   
      monomers[n1].force_areaL[j]=0.; 
      monomers[n1].force_areaG[j]=0.;
      monomers[n1].force_wall[j]=0.;  
      monomers[n1].force_inter[j]=0.;
      monomers[n1].force_interA[j]=0.;
      monomers[n1].force_interR[j]=0.;              
		}
	}
	for(int n1=0; n1 < numBead; n1++) {
		for(int i=0; i < DIMS; i++) {
			for(int j=0; j < DIMS; j++) {
				monomers[n1].stress[i][j]=0.;
        monomers[n1].stress_elas[i][j]=0.;
        monomers[n1].stress_bend[i][j]=0.;
        monomers[n1].stress_vol[i][j]=0.;
        monomers[n1].stress_areaL[i][j]=0.;
        monomers[n1].stress_areaG[i][j]=0.; 
        monomers[n1].stress_wall[i][j]=0.;
        monomers[n1].stress_int_v1[i][j]=0.;
        monomers[n1].stress_int_v2[i][j]=0.;
			}
		}
	}
}  

void get_forces(struct sphere_param *sphere_pm, struct monomer *monomers, struct face 
                *faces, Float ***velcs_df, int n_step, VSLStreamStatePtr rngstream)
{
	extern int n_proc, num_proc;
	extern int num_x;
	extern int max_x, max_y, max_z;
	extern int wall_flag;
	int mon_proc;
	int num_beads;
	int i, j;
	int n1, n2;
	int s_type = 0;              /* spring type : 0=FENE, 1=WLC, 2=harmonic */
	int ev_type = 1;             /* ev type: 0=HS, 1=WCA, 2=gaussian */
	int maxsize[DIMS];
	int max_N_flag=0;
	double send_buffer[2*DIMS+1];
	char filename[80];
	FILE *stream;
	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;
	num_beads = sphere_pm->num_beads;
  extern char *work_dir;

	//for(n1=0; n1 < num_beads; n1++) 
  //{
  //  monomers[n1].updatedFlag=FALSE; // Modification 20170725  
	//	for(j=0; j < DIMS; j++) { 
	//		monomers[n1].force[j] = 0.0;
  //    monomers[n1].force_ev[j]=0.;
  //    monomers[n1].force_face[j]=0.;
  //    monomers[n1].force_spring[j]=0.;
  //    monomers[n1].force_bend[j]=0.;
  //    monomers[n1].force_vol[j]=0.;   
  //    monomers[n1].force_area[j]=0.; 
  //    monomers[n1].drag[j]=0.;      
  //    monomers[n1].force_wall[j]=0.;  
  //    monomers[n1].force_fluc[j]=0.;        
  //    monomers[n1].force_pp[j] = 0.0;
  //    monomers[n1].force_nonbonded_intra[j]=0.;              
	//	}
	//}
	//for(n1=0; n1 < num_beads; n1++) {
	//	for(i=0; i < DIMS; i++) {
	//		for(j=0; j < DIMS; j++) {
	//			monomers[n1].stress[i][j] = 0.0;
  //      monomers[n1].stress_elas[i][j] = 0.0;
  //      monomers[n1].stress_bend[i][j] = 0.0;
  //      monomers[n1].stress_vol[i][j] = 0.0;
  //      monomers[n1].stress_area[i][j] = 0.0;
  //      monomers[n1].stress_wall[i][j] = 0.0;
  //      monomers[n1].stress_int_v1[i][j] = 0.0;
  //      monomers[n1].stress_int_v2[i][j] = 0.0;
	//		}
	//	}
	//}

  zero_variable(sphere_pm->num_beads, monomers);

  /* Modification 20170316: write force information
     --------------------- 
  */
//  sprintf(filename,"%s/data/forceInfor.dat",work_dir); // TEMPORARY!
//  stream=fopen(filename,"w");
//  fprintf(stream, "step = %d\n", n_step);
//  fclose(stream);

	// Calculate spring and exc. vol. forces if # of beads > 1
	if(num_beads > 1) 
  {
		spring_force(sphere_pm, monomers, faces);
 
    // Modification 20170323
		if(num_beads < 200) {
      //nonbonded_interaction(sphere_pm, monomers);  // 20170308
	    ev_force(monomers[0].radius, sphere_pm, monomers);	
    }	
		else 
    {
			for(n1=0; n1 < num_beads; n1++) {
				if(monomers[n1].list[0] >= MAX_N-1) 
					max_N_flag = 1;				
			}
			if(max_N_flag == 1) 
		    ev_force(monomers[0].radius, sphere_pm, monomers);			
			else 
  	    ev_force_nlist(monomers[0].radius, sphere_pm, monomers);
    }
	}
	for(n1=0; n1 < num_beads; n1++) {
		for(j=0; j<DIMS; j++) {
			//monomers[n1].force[j] *= sphere_pm->monmass; // why?? Modification 20170419
			monomers[n1].f_int[j] = monomers[n1].force[j];
		}
	}

	if(wall_flag > 0) 
    wall_force(sphere_pm, monomers); 

	/* send forces and fluid_vel back to proc 0 then broadcast from 0 */
	for(n1=0; n1 < num_beads; n1++) { 
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

	// Calculate particle stress. The contribution from inter-particle
	// interaction is calculated in 'ev_force' function.
	//ParticleStress (sphere_pm, monomers);
  for(int n=0; n < num_beads; n++)
    for(i=0; i < DIMS; i++)
      for(j=0; j < DIMS; j++)
        monomers[n].stress[i][j] += ( monomers[n].stress_elas[i][j] + 
        monomers[n].stress_bend[i][j] + monomers[n].stress_vol[i][j] + 
        monomers[n].stress_areaG[i][j] + monomers[n].stress_wall[i][j] );
}

void get_force_growth(struct sphere_param *sphere_pm, struct monomer *monomers, struct face 
     *faces, Float ***velcs_df, int n_step, VSLStreamStatePtr rngstream)
{
	extern int n_proc, num_proc, num_x, max_x, max_y, max_z, wall_flag;
	int mon_proc;
	int num_beads;
	int i, j;
	int n1, n2;
	int s_type = 0;              /* spring type : 0=FENE, 1=WLC, 2=harmonic */
	int ev_type = 1;             /* ev type: 0=HS, 1=WCA, 2=gaussian */
	int maxsize[DIMS];
	int max_N_flag=0;
	double send_buffer[2*DIMS+1];
	char filename[80];
	FILE *stream;
	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;
	num_beads = sphere_pm->num_beads;
  extern char *work_dir;

	//for(n1=0; n1 < num_beads; n1++) 
  //{
	//	for(j=0; j < DIMS; j++) { 
	//		monomers[n1].force[j] = 0.0;
  //    monomers[n1].force_ev[j]=0.;
  //    monomers[n1].force_face[j]=0.;
  //    monomers[n1].force_spring[j]=0.;
  //    monomers[n1].force_bend[j]=0.;
  //    monomers[n1].force_vol[j]=0.;   
  //    monomers[n1].force_areaG[j]=0.; 
  //    monomers[n1].drag[j]=0.;      
  //    monomers[n1].force_wall[j]=0.;  
  //    monomers[n1].force_fluc[j]=0.;        
  //    monomers[n1].force_pp[j] = 0.0;
  //    monomers[n1].force_nonbonded_intra[j]=0.;              
	//	}
	//}
	//for(n1=0; n1 < num_beads; n1++) {
	//	for(i=0; i < DIMS; i++) {
	//		for(j=0; j < DIMS; j++) {
	//			monomers[n1].stress[i][j] = 0.0;
  //      monomers[n1].stress_elas[i][j] = 0.0;
  //      monomers[n1].stress_bend[i][j] = 0.0;
  //      monomers[n1].stress_vol[i][j] = 0.0;
  //      monomers[n1].stress_areaG[i][j] = 0.0;
  //      monomers[n1].stress_wall[i][j] = 0.0;
  //      monomers[n1].stress_int_v1[i][j] = 0.0;
  //      monomers[n1].stress_int_v2[i][j] = 0.0;
	//		}
	//	}
	//}
  zero_variable(sphere_pm->num_beads,monomers);

  /* Modification 20170316: write force information
     --------------------- 
  */
//  sprintf(filename,"%s/data/forceInfor.dat",work_dir); // TEMPORARY!
//  stream=fopen(filename,"w");
//  fprintf(stream, "step = %d\n", n_step);
//  fclose(stream);

	// Calculate spring and exc. vol. forces if # of beads > 1
	if(num_beads > 1) 
  {
		spring_force(sphere_pm, monomers, faces);
 
    // Modification 20170323
		if(num_beads < 200) {
      //nonbonded_interaction(sphere_pm, monomers);  // 20170308
	    ev_force(monomers[0].radius, sphere_pm, monomers);	
    }	
		else 
    {
			for(n1=0; n1 < num_beads; n1++) {
				if(monomers[n1].list[0] >= MAX_N-1) 
					max_N_flag = 1;				
			}
			if(max_N_flag == 1) 
		    ev_force(monomers[0].radius, sphere_pm, monomers);			
			else 
  	    ev_force_nlist(monomers[0].radius, sphere_pm, monomers);
    }
    
    nonbonded_interaction_nlist(sphere_pm, monomers); // Modification 20170319
	}

	if(wall_flag > 0) 
    wall_force(sphere_pm, monomers); 

	/* send forces and fluid_vel back to proc 0 then broadcast from 0 */
	for(n1=0; n1 < num_beads; n1++) { 
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

void get_forces_hi(struct sphere_param *sphere_pm, struct monomer *monomers, struct face 
                *faces, Float ***velcs_df, int n_step, VSLStreamStatePtr rngstream)
{
	extern int n_proc, num_proc;
	extern int num_x;
	extern int max_x, max_y, max_z;
	extern int wall_flag;
	int mon_proc;
	int num_beads;
	int i, j;
	int n1, n2;
	int s_type = 0;              /* spring type : 0=FENE, 1=WLC, 2=harmonic */
	int ev_type = 1;             /* ev type: 0=HS, 1=WCA, 2=gaussian */
	int maxsize[DIMS];
	int max_N_flag=0;
	double send_buffer[2*DIMS+1];
	char filename[80];
	FILE *stream;
	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;
	num_beads = sphere_pm->num_beads;
  extern char *work_dir;

	//  for(j=0; j<DIMS; j++) 
	//    monomers[0].force[j] = monomers[0].f_ext[j];

	for(n1=0; n1 < num_beads; n1++) {
		for(j=0; j < DIMS; j++) { 
			monomers[n1].force[j] = 0.0;
      monomers[n1].force_ev[j]=0.;
      monomers[n1].force_face[j]=0.;
      monomers[n1].force_spring[j]=0.;
      monomers[n1].force_bend[j]=0.;
      monomers[n1].force_vol[j]=0.;   
      monomers[n1].force_areaG[j]=0.; 
      monomers[n1].drag[j]=0.;      
      monomers[n1].force_wall[j]=0.;  
      monomers[n1].force_fluc[j]=0.;        
      monomers[n1].force_pp[j] = 0.0;
      monomers[n1].force_nonbonded_intra[j]=0.;              
		}
	}
	for(n1=0; n1 < num_beads; n1++) {
		for(i=0; i < DIMS; i++) {
			for(j=0; j < DIMS; j++) {
				monomers[n1].stress[i][j] = 0.0;
        monomers[n1].stress_elas[i][j] = 0.0;
        monomers[n1].stress_bend[i][j] = 0.0;
        monomers[n1].stress_vol[i][j] = 0.0;
        monomers[n1].stress_areaG[i][j] = 0.0;
        monomers[n1].stress_wall[i][j] = 0.0;
        monomers[n1].stress_int_v1[i][j] = 0.0;
        monomers[n1].stress_int_v2[i][j] = 0.0;
			}
		}
	}

  /* Modification 20170316: write force information
     --------------------- 
  */
//  sprintf(filename,"%s/data/forceInfor.dat",work_dir); // TEMPORARY!
//  stream=fopen(filename,"w");
//  fprintf(stream, "step = %d\n", n_step);
//  fclose(stream);

	// Calculate spring and exc. vol. forces if # of beads > 1
	if(num_beads > 1) 
  {
		spring_force(sphere_pm, monomers, faces);
 
    // Modification 20170323
		if(num_beads < 200) {
      //nonbonded_interaction(sphere_pm, monomers);  // 20170308
	    ev_force(monomers[0].radius, sphere_pm, monomers);	
    }	
		else 
    {
			for(n1=0; n1 < num_beads; n1++) {
				if(monomers[n1].list[0] >= MAX_N-1) 
					max_N_flag = 1;				
			}
			if(max_N_flag == 1) 
		    ev_force(monomers[0].radius, sphere_pm, monomers);			
			else 
  	    ev_force_nlist(monomers[0].radius, sphere_pm, monomers);
        //ev_force(monomers[0].radius, sphere_pm, monomers);  //test
    }
    
    nonbonded_interaction_nlist(sphere_pm, monomers); // Modification 20170319
    /*if(n_step > 2000) // Note: should be modified. 20161003
      AggreForce2(sphere_pm, monomers, work_dir);
    */ //20161129
	}
	for(n1=0; n1 < num_beads; n1++) {
		for(j=0; j<DIMS; j++) {
			//monomers[n1].force[j] *= sphere_pm->monmass; // why?? Modification 20170419
			monomers[n1].f_int[j] = monomers[n1].force[j];
		}
	}
  // Modification 20170419
	hi_force(monomers, velcs_df, sphere_pm->fric, sphere_pm->MD_steps, sphere_pm->dt, 
           num_beads, sphere_pm->monmass, rngstream);

	if(wall_flag > 0) 
    wall_force(sphere_pm, monomers); 

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

	// Calculate particle stress. The contribution from inter-particle
	// interaction is calculated in 'ev_force' function.
	//ParticleStress (sphere_pm, monomers);
  for(int n=0; n < num_beads; n++)
    for(i=0; i < DIMS; i++)
      for(j=0; j < DIMS; j++)
        monomers[n].stress[i][j] += ( monomers[n].stress_elas[i][j] + 
        monomers[n].stress_bend[i][j] + monomers[n].stress_vol[i][j] + 
        monomers[n].stress_areaG[i][j] + monomers[n].stress_wall[i][j] );
}


/* Calculate spring force : type = 0 FENE, 1 WLC, 2 harmonic */
/* Calculate volume constraint force */
void spring_force(struct sphere_param *sphere_pm, struct monomer *monomers, struct face 
     *faces)
{
	extern int n_proc, max_x, max_y, max_z, num_x, wall_flag;
	int n1, n2, n3, n4;
	int i,j,k,d;
	int num_beads = sphere_pm->num_beads;
	int Nsphere[NTYPES];
	int Nface;
	int Nbead;
	int nlevel[NTYPES];
	int maxsize[DIMS];
	int nface, numsphere, type, face0, bead0;
	Float r1;
	Float slope, intercept;
	Float mag;
	Float Nks = sphere_pm->nks*sphere_pm->sigma_k;
	Float sigma = monomers[0].radius*2.0;
	Float q12[DIMS], q1[DIMS], q2[DIMS], q3[DIMS];
	Float q0, q12mag;
	Float force[DIMS];
	double normal[2][DIMS], mag_normal[2], a[DIMS];
	double r_face_c[DIMS], dr_face_c[3][DIMS], mag_dr_c[3];
	double dot1a, dot2a, dot3a;
	double n12, n11, n22;
	double H_fene[NTYPES];
	double k_bend[NTYPES];
	double k_V[NTYPES], k_A[NTYPES];
	double dr[DIMS];
	double V0[NTYPES], A0[NTYPES];
	double **com, *delta_V, *delta_A;

  extern char *work_dir;
  char fileName[100];
	FILE *stream;
  
  extern double springLength_eq;  // modified

	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;
	for(j=0; j < NTYPES; j++) {
		Nsphere[j] = sphere_pm->Ntype[j];
		nlevel[j] = sphere_pm->nlevel[j];
		H_fene[j] = sphere_pm->H_fene[j];
		k_bend[j] = sphere_pm->k_bend[j]*sphere_pm->kT/sigma;
		k_V[j] = sphere_pm->k_V[j]*sphere_pm->kT/sigma;
		k_A[j] = sphere_pm->k_A[j]*sphere_pm->kT/sigma;
		V0[j] = sphere_pm->V0[j];
		A0[j] = sphere_pm->A0[j];
	}
	numsphere = Nsphere[0] + Nsphere[1];

	com = (double **) calloc(numsphere, sizeof(double *));
	delta_V = (double *) calloc(numsphere, sizeof(double));
	delta_A = (double *) calloc(numsphere, sizeof(double));
	for(i=0; i<numsphere; i++)
		com[i] = (double *) calloc(DIMS, sizeof(double));

  // Loop over all particles
	for(i=0; i < numsphere; i++) 
  {
		type = (i<sphere_pm->Ntype[0] ? 0 :1);
		Nbead = sphere_pm->N_per_sphere[type];
		Nface = sphere_pm->face_per_sphere[type];
		face0 = (type==0 ? i*sphere_pm->face_per_sphere[0] : 
            Nsphere[0]*sphere_pm->face_per_sphere[0] + 
            (i-Nsphere[0])*sphere_pm->face_per_sphere[1]);
		bead0 = (type==0 ? i*sphere_pm->N_per_sphere[0] : 
            Nsphere[0]*sphere_pm->N_per_sphere[0] + 
            (i-Nsphere[0])*sphere_pm->N_per_sphere[1]);

		// Calculate the center-of-mass of a sphere
		for(d=0; d < DIMS; d++) {
			com[i][d] = 0.0;
			for(j=0; j < Nbead; j++)
				com[i][d] += monomers[bead0+j].pos_pbc[d];
			com[i][d] /= Nbead;
		}
		/* calculate face area and total volume */
		/* loop over all faces */
		delta_V[i] =0.;
		delta_A[i] =0.;
		for(j=0; j < sphere_pm->face_per_sphere[type]; j++) 
    {
			nface = face0 + j;
			n1 = faces[nface].vertices[0];
			for(d=0; d < DIMS; d++)
				dr[d] = monomers[n1].pos_pbc[d] - com[i][d];
			n2 = faces[nface].vertices[1];
			for(d=0; d < DIMS; d++)
				q1[d] = monomers[n2].pos_pbc[d] - monomers[n1].pos_pbc[d];
			n3 = faces[nface].vertices[2];
			for(d=0; d < DIMS; d++)
				q2[d] = monomers[n3].pos_pbc[d] - monomers[n1].pos_pbc[d];
			product(q1, q2, normal[0]);
			delta_V[i] += fabs(dr[0]*normal[0][0]+dr[1]*normal[0][1]+dr[2]*normal[0][2])/6.;
			delta_A[i] += sqrt(normal[0][0]*normal[0][0]+normal[0][1]*normal[0][1]+
                         normal[0][2]*normal[0][2])/2.0;
		}
		/* substract the equilibrium volume */
		delta_V[i] = (-1)*k_V[type]*(delta_V[i]-V0[type])/V0[type];
		delta_A[i] = (-1)*k_A[type]*(delta_A[i]-A0[type])/A0[type];
		/* force is divided upon 3 nodes. another factor of 2 due to area */
		delta_V[i] /= 6.;
    
    // Compute elastic and bending forces; Loop over all edges 
		for(j=0; j < Nbead; j++) 
    {
			n1 = bead0 + j;

			for(k=1; k <= monomers[n1].blist[0][0] ; k++) 
      {
				n2 = monomers[n1].blist[k][0];

				q12mag=0.0;
				for(d=0; d < DIMS; d++) {
				  q12[d] = monomers[n1].pos_pbc[d] - monomers[n2].pos_pbc[d];
					//if((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1)) // Modification 20170811: unnecessary.
					//  q12[d]=n_image(q12[d], maxsize[d]);
					q12mag += q12[d]*q12[d];
				}    
				q12mag = sqrt(q12mag);
				q12[0] /= q12mag;    q12[1] /= q12mag;    q12[2] /= q12mag;

				if(sphere_pm->spring == 0) // FENE spring
        {
					q0 = sphere_pm->Q_fene[type]*(monomers[n1].radius+monomers[n2].radius);
					if(q12mag < q0) 
						mag = (sphere_pm->kT/sigma)*H_fene[type]*q12mag/(1.0 - (q12mag/q0)*(q12mag/q0));
					else 
          {
						stream = fopen("output.err", "w");
						fprintf(stream, "Spring is broken, %d %d %le %le\n", n1, n2, q12mag, q0);
						fprintf(stream, "sphere %d mon %d pos=(%le %le %le), mon %d pos = (%le %le %le)\n", 
                i, n1, monomers[n1].pos_pbc[0], monomers[n1].pos_pbc[1], 
                monomers[n1].pos_pbc[2], n2, monomers[n2].pos_pbc[0], 
                monomers[n2].pos_pbc[1], monomers[n2].pos_pbc[2]);
						fprintf(stream, "mon %d force=(%le %le %le), mon %d force = (%le %le %le)\n", 
                n1, monomers[n1].force[0], monomers[n1].force[1], monomers[n1].force[2],
                n2, monomers[n2].force[0], monomers[n2].force[1], monomers[n2].force[2]);
						for(k=1; k < monomers[n1].blist[0][0]; k++) {
							n2 = monomers[n1].blist[k][0];
							fprintf(stream, "%d %le %le %le\n", n2, monomers[n2].pos_pbc[0], 
                  monomers[n2].pos_pbc[1], monomers[n2].pos_pbc[2]);
						}
						fprintf(stream, "\n");
						fclose(stream);
						exit(11);
					}
				}
				else if(sphere_pm->spring == 1) // Worm like spring
        {
					//r1 = q12mag / Nks; // modified 1020
          double springLength_max;
          r1 = q12mag / springLength_max;
					if(r1 < 1.0) {
            double temp = 1.0 - r1;  // modified 1020
						//mag = 2.0*sphere_pm->kT/sphere_pm->sigma_k*(0.25/(temp*temp)- 0.25 + r1);
            mag = sphere_pm->kT / sphere_pm->persisLength * (0.25/(temp*temp)- 0.25 + r1);
          }
					else {	  
						if(n_proc == 0) {
							stream = fopen("output.err", "w");
							fprintf(stream, "Spring is broken, %d %le\n", n1, q12mag);
							fprintf(stream, "sphere %d mon %d pos=(%le %le %le), mon %d pos = (%le %le %le)\n", 
                  i, n1, monomers[n1].pos_pbc[0], monomers[n1].pos_pbc[1], 
                  monomers[n1].pos_pbc[2], n2, monomers[n2].pos_pbc[0], 
                  monomers[n2].pos_pbc[1], monomers[n2].pos_pbc[2]);
							fprintf(stream, "mon %d force=(%le %le %le), mon %d force = (%le %le %le)\n", 
                  n1, monomers[n1].force[0], monomers[n1].force[1], monomers[n1].force[2], 
                  n2, monomers[n2].force[0], monomers[n2].force[1], monomers[n2].force[2]);
							fclose(stream);
							exit(11);
						}
					}
				}
				else if(sphere_pm->spring == 2) // harmonic spring
        {
          // Note: mag = (partial U / partial r) 
  				//mag = 2.0*(sphere_pm->kT/sigma)*H_fene[type]*(q12mag-1.0);
          mag = 2.0*sphere_pm->kT*H_fene[type]*(q12mag-monomers[n1].initLength[k]);
        }

				for(d=0; d < DIMS; d++) {
					force[d]= q12[d] * mag;            // assign the direction of the force
					monomers[n1].force[d] -= force[d]; // f = -(partial U / partial r)
					monomers[n2].force[d] += force[d];
          monomers[n1].force_spring[d] -= force[d];
          monomers[n2].force_spring[d] += force[d];
        }
        for(int m=0; m < DIMS; m++) {
          for(int n=0; n < DIMS; n++) {
            monomers[n1].stress_elas[m][n] -= monomers[n1].pos[m]*force[n];
            monomers[n2].stress_elas[m][n] += monomers[n2].pos[m]*force[n];
          }
        }

        // Bending force on face (n1,n2,n3) and face (n1,n2,n4)    20161124
        double x12[DIMS], x21[DIMS], x13[DIMS], x31[DIMS], x24[DIMS], x32[DIMS], x41[DIMS];
        double normal1[DIMS], normal2[DIMS];
        double normal1_mag, normal2_mag;
        double cos_theta, sin_theta;
        double normal1_minus_cos_theta_normal2[DIMS];
        double normal2_minus_cos_theta_normal1[DIMS];
        double factor;
        double term11[DIMS], term12[DIMS], term21[DIMS], term22[DIMS], term3[DIMS], term4[DIMS];
        double f1[DIMS], f2[DIMS], f3[DIMS], f4[DIMS];
        double kBend = 2./sqrt(3)*k_bend[0];//sphere_pm->kBend; // convert k_c to k_b. temporary setting
        double theta0; 
        //if(nlevel[0]==2 || nlevel[0]==-1)   
        //  //theta0 = 0.2171;//sphere_pm->theta0; only temporary setting
        //  theta0 = monomers[n1].initAngle[k]; //20170219
        //if(nlevel[0]==3 || nlevel[0]==-1) 
        //  //theta0 = 0.1070;  // 20170202
        theta0 = monomers[n1].initAngle[k]; //20170219  // Modification 20170809

  
        double x34[DIMS], normal_12[DIMS]; // 20170124
 
        n3 = monomers[n1].blist[k][1]; // n1 = bead0+j; n2 = monomers[n1].blist[k][0];
        n4 = monomers[n1].blist[k][2];
        for(d=0; d < DIMS; d++) {
          x12[d] = monomers[n1].pos_pbc[d] - monomers[n2].pos_pbc[d];
          x21[d] = -x12[d];
          x13[d] = monomers[n1].pos_pbc[d] - monomers[n3].pos_pbc[d];
          x31[d] = -x13[d];
          x24[d] = monomers[n2].pos_pbc[d] - monomers[n4].pos_pbc[d];
          x32[d] = monomers[n3].pos_pbc[d] - monomers[n2].pos_pbc[d];
          x41[d] = monomers[n4].pos_pbc[d] - monomers[n1].pos_pbc[d];
        }
        product(x31, x21, normal1);
        product(x21, x41, normal2);
        normal1_mag = sqrt(normal1[0]*normal1[0] + normal1[1]*normal1[1] + 
                      normal1[2]*normal1[2]);
        normal2_mag = sqrt(normal2[0]*normal2[0] + normal2[1]*normal2[1] +
                      normal2[2]*normal2[2]);
        normal1[0] /= normal1_mag;  normal1[1] /= normal1_mag;  normal1[2] /= normal1_mag;
        normal2[0] /= normal2_mag;  normal2[1] /= normal2_mag;  normal2[2] /= normal2_mag;
    
        //if(isnan(normal1[0])) 
        //{
        //  sprintf(fileName,"%s/data/nan_warning.dat",work_dir);
        //  stream=fopen(fileName,"a");
        //  fprintf(stream,"normal1 =(%f, %f, %f) normal1_mag=%f\n",normal1[0],normal1[1],
        //          normal1[2],normal1_mag);
        //  fclose(stream);
        //  exit(18);          
        //}
        //if(isnan(normal2[0])) 
        //{
        //  sprintf(fileName,"%s/data/nan_warning.dat",work_dir);
        //  stream=fopen(fileName,"a");
        //  fprintf(stream,"normal2 =(%f, %f, %f) normal2_mag=%f\n",normal2[0],normal2[1],
        //          normal2[2],normal2_mag);
        //  fclose(stream);
        //  exit(18);          
        //}

        double n1_cross_n2[3];
        product(normal1,normal2,n1_cross_n2);
        double norm_n1crossn2 = sqrt(n1_cross_n2[0]*n1_cross_n2[0]+
               n1_cross_n2[1]*n1_cross_n2[1] + n1_cross_n2[2]*n1_cross_n2[2]);
        double theta;          
        theta = atan2(norm_n1crossn2, iproduct(normal1,normal2));

        //sprintf(fileName,"%s/data/angle.dat",work_dir);
        //stream=fopen(fileName,"a");
        //fprintf(stream, "angle = %f\n", theta*180/Pi);
        //fclose(stream);

//        cos_theta = iproduct(normal1, normal2);  Modification 20170322
        cos_theta = cos(theta); 
               
        // version 2 determine sin(theta)
/*        for(d=0; d<DIMS; d++) {
          normal_12[d] = normal1[d] - normal2[d];
          x34[d] = monomers[n3].pos_pbc[d] - monomers[n4].pos_pbc[d];
        }
        if(iproduct(normal_12, x34) >=0.) {
          sin_theta = sqrt(1. - cos_theta * cos_theta);
        }
        else
          sin_theta = -sqrt(1. - cos_theta * cos_theta);
*/

          //Modification 20170322  
//        sin_theta = sqrt(1. - cos_theta * cos_theta);  // version 1
        sin_theta = sin(theta); 

        //if(isnan(sin_theta)) 
        //{
        //  sprintf(fileName,"%s/data/nan_warning.dat",work_dir);
        //  stream=fopen(fileName,"a");
        //  fprintf(stream,"cos_theta = %le\n",cos_theta);
        //  fclose(stream);
        //  exit(18);          
        //}

        for(d=0; d < DIMS; d++) {
          normal1_minus_cos_theta_normal2[d] = normal1[d] - cos_theta*normal2[d];
          normal2_minus_cos_theta_normal1[d] = normal2[d] - cos_theta*normal1[d];
        }
        factor = kBend * (sin_theta*cos(theta0)-cos_theta*sin(theta0)) / sin_theta;

        //if(isnan(factor)) 
        //{
        //  sprintf(fileName,"%s/data/nan_warning.dat",work_dir);
        //  stream=fopen(fileName,"a");
        //  fprintf(stream,"sin_theta = %f\n",sin_theta);
        //  fclose(stream);
        //  exit(18);          
        //}

        product(x24, normal1_minus_cos_theta_normal2, term11);
        product(x32, normal2_minus_cos_theta_normal1, term12);
        product(x41, normal1_minus_cos_theta_normal2, term21);
        product(x13, normal2_minus_cos_theta_normal1, term22);
        product(x21, normal2_minus_cos_theta_normal1, term3);
        product(x12, normal1_minus_cos_theta_normal2, term4);
        for(d=0; d < DIMS; d++) {
          f1[d] = factor*(term11[d] / normal2_mag + term12[d] / normal1_mag);
          f2[d] = factor*(term21[d] / normal2_mag + term22[d] / normal1_mag);
          f3[d] = factor*term3[d] / normal1_mag;
          f4[d] = factor*term4[d] / normal2_mag;
        }
        for(d=0 ; d<DIMS ; d++) {
          monomers[n1].force[d] += f1[d];
          monomers[n2].force[d] += f2[d];
          monomers[n3].force[d] += f3[d];
          monomers[n4].force[d] += f4[d];
          monomers[n3].force_bend[d] += f1[d];
          monomers[n4].force_bend[d] += f2[d];
          monomers[n1].force_bend[d] += f3[d];
          monomers[n2].force_bend[d] += f4[d];
        }
//        /* Modification 20170316: write force information
//           ---------------------
//        */
//        double f1mag = sqrt(f1[0]*f1[0]+f1[1]*f1[1]+f1[2]*f1[2]);
//        double f2mag = sqrt(f2[0]*f2[0]+f2[1]*f2[1]+f2[2]*f2[2]);
//        double f3mag = sqrt(f3[0]*f3[0]+f3[1]*f3[1]+f3[2]*f3[2]);
//        double f4mag = sqrt(f4[0]*f4[0]+f4[1]*f4[1]+f4[2]*f4[2]);
//        if(f1mag > 0.3 || f2mag > 0.3 || f3mag > 0.3 ||f4mag > 0.3) 
//        {
//          sprintf(fileName,"%s/data/warningMes.dat",work_dir); 
//          stream=fopen(fileName,"a");
//          fprintf(stream, "bending force = (%f, %f, %f, %f)\n", f1mag, f2mag, f3mag, f4mag);
//          fclose(stream);
//        }
        for(int m=0; m < DIMS; m++) {
          for(int n=0; n < DIMS; n++) {
            monomers[n1].stress_bend[m][n] += monomers[n1].pos[m]*f1[n];
            monomers[n2].stress_bend[m][n] += monomers[n2].pos[m]*f2[n];
            monomers[n3].stress_bend[m][n] += monomers[n3].pos[m]*f3[n];
            monomers[n4].stress_bend[m][n] += monomers[n4].pos[m]*f4[n];
          }
        }   
			}
		}
		// Compute global and local area restoring forces and volume constraint force
		for(j=0; j < Nface; j++) 
    {
			nface = face0 + j;
			n1 = faces[nface].vertices[0];
			n2 = faces[nface].vertices[1];
			n3 = faces[nface].vertices[2];
      double v23[DIMS], v31[DIMS], v12[DIMS], v21[DIMS];
      double normalVector[DIMS], normal_mag, normalizedNormal[DIMS];
      double temp_n1[DIMS], temp_n2[DIMS], temp_n3[DIMS];
      double faceCenter[DIMS];            
      double temp_v1[DIMS], temp_v2[DIMS], temp_v3[DIMS]; 

      for(d=0; d < DIMS; d++) 
        faceCenter[d] = (monomers[n1].pos_pbc[d] + monomers[n2].pos_pbc[d] + 
                         monomers[n3].pos_pbc[d]) / 3.;
      for(d=0; d < DIMS; d++) {
        v23[d] = monomers[n2].pos_pbc[d] - monomers[n3].pos_pbc[d];
        v31[d] = monomers[n3].pos_pbc[d] - monomers[n1].pos_pbc[d];
        v12[d] = monomers[n1].pos_pbc[d] - monomers[n2].pos_pbc[d];
        v21[d] = monomers[n2].pos_pbc[d] - monomers[n1].pos_pbc[d];
      }
			product(v31, v21, normalVector);
      normalizedNormal[0] = normalVector[0];
      normalizedNormal[1] = normalVector[1];
      normalizedNormal[2] = normalVector[2];
      normal_mag = sqrt(iproduct(normalVector, normalVector));
      normalizedNormal[0] /= normal_mag;
      normalizedNormal[1] /= normal_mag;
      normalizedNormal[2] /= normal_mag;
      double triangleArea = 0.5 * normal_mag;			
			product(normalizedNormal, v23, temp_n1);
      product(normalizedNormal, v31, temp_n2);
      product(normalizedNormal, v12, temp_n3);
      product(faceCenter, v23, temp_v1);
      product(faceCenter, v31, temp_v2);
      product(faceCenter, v12, temp_v3); 
      // Global area constraint
      double force1[DIMS], force2[DIMS], force3[DIMS];
      for(d=0; d<DIMS; d++)
      {
        force1[d] = delta_A[i]*0.5*temp_n1[d];
        force2[d] = delta_A[i]*0.5*temp_n2[d];
        force3[d] = delta_A[i]*0.5*temp_n3[d];
      }
			for(d=0 ; d < DIMS ; d++) {
		    monomers[n1].force[d] += force1[d];
				monomers[n2].force[d] += force2[d];
				monomers[n3].force[d] += force3[d];
        monomers[n1].force_areaG[d] += force1[d];
        monomers[n2].force_areaG[d] += force2[d];
        monomers[n3].force_areaG[d] += force3[d];
      }
      for(int m=0; m < DIMS; m++) {
        for(int n=0; n < DIMS; n++) {
          monomers[n1].stress_areaG[m][n] += monomers[n1].pos[m] * force1[n];
          monomers[n2].stress_areaG[m][n] += monomers[n2].pos[m] * force2[n];
          monomers[n3].stress_areaG[m][n] += monomers[n3].pos[m] * force3[n];
		    }
	    } 
      // Local area constraint
      double factor = -0.25 * sphere_pm->ka_local * (triangleArea - faces[nface].area_0) /
                      (triangleArea * faces[nface].area_0);
      product(normalVector, v23, temp_n1);
      product(normalVector, v31, temp_n2);
      product(normalVector, v12, temp_n3);
      force1[0] = factor * temp_n1[0];
      force2[0] = factor * temp_n2[0];
      force3[0] = factor * temp_n3[0];
      force1[1] = factor * temp_n1[1];
      force2[1] = factor * temp_n2[1];
      force3[1] = factor * temp_n3[1];
      force1[2] = factor * temp_n1[2];
      force2[2] = factor * temp_n2[2];
      force3[2] = factor * temp_n3[2];
      for(d=0; d<DIMS; d++) {
        monomers[n1].force[d] += force1[d];
        monomers[n2].force[d] += force2[d];
        monomers[n3].force[d] += force3[d];
        monomers[n1].force_areaL[d] += force1[d];
        monomers[n2].force_areaL[d] += force2[d];
        monomers[n3].force_areaL[d] += force3[d];
      }
      for(int m=0; m < DIMS; m++) {
        for(int n=0; n < DIMS; n++) {
          monomers[n1].stress_areaL[m][n] += monomers[n1].pos[m] * force1[n];
          monomers[n2].stress_areaL[m][n] += monomers[n2].pos[m] * force2[n];
          monomers[n3].stress_areaL[m][n] += monomers[n3].pos[m] * force3[n];
		    }
	    } 
      // Volume constraint
      for(d=0; d<DIMS; d++)
      {
        force1[d] = delta_V[i]*(normalVector[d]/3. + temp_v1[d]);
        force2[d] = delta_V[i]*(normalVector[d]/3. + temp_v2[d]);
        force3[d] = delta_V[i]*(normalVector[d]/3. + temp_v3[d]);
      }
      for(d=0; d<DIMS; d++) {
 		    monomers[n1].force[d] += force1[d];
				monomers[n2].force[d] += force2[d];
				monomers[n3].force[d] += force3[d];
        monomers[n1].force_vol[d] += force1[d];
        monomers[n2].force_vol[d] += force2[d];
        monomers[n3].force_vol[d] += force3[d];       
			}
      for(int m=0; m < DIMS; m++) {
        for(int n=0; n < DIMS; n++) {
          monomers[n1].stress_vol[m][n] += monomers[n1].pos[m] * force1[n];
          monomers[n2].stress_vol[m][n] += monomers[n2].pos[m] * force2[n];
          monomers[n3].stress_vol[m][n] += monomers[n3].pos[m] * force3[n];
		    }
	    } 
    }
  }
	free(delta_V);
	free(delta_A);
	for(i=0; i<numsphere; i++)
		free(com[i]);
	free(com); 
}

/* Exc. Vol forces type = 0 HS, 1 LJ, 2 Gaussian */
void ev_force(Float radius, struct sphere_param *sphere_pm, struct monomer *monomers)
{
	extern int n_proc;
	extern int max_x, max_y, max_z;
	extern int num_x;
	extern int wall_flag;
	extern double tau;
	int i, j;
	long n1, n2;
	int num_beads = sphere_pm->num_beads;
	Float sigma = 2.0*radius;
	Float cutoff = sphere_pm->evcutoff*radius;
	Float aterm;
	Float q12[DIMS];
	Float q12mag;
	Float eps, Ecut;
	Float force[DIMS];
	int maxsize[DIMS];
	FILE *stream;
	double eta = (tau - 0.5)/3.0;
	double Rg = sphere_pm->Ss;
	double Nks =sphere_pm->nks;
	double v = sphere_pm->sigma_k*sphere_pm->sigma_k*sphere_pm->sigma_k;
	double aev_alpha = 4.0/3.0*Rg*Rg;
	double aev = Nks*Nks*(1.0/((Pi*aev_alpha)*sqrt(Pi*aev_alpha)));
	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;

  extern double springLength_eq;

  sigma = 1 / 1.122;
	cutoff = 1.0;
  //sigma = springLength_eq / 1.122;
  //cutoff = springLength_eq;
//	eps = sphere_pm->kT/sigma;
  eps = sphere_pm->epsilon_WCA * sphere_pm->kT;

	for(i=0; i<DIMS; i++)
		force[i]=0.0;

	for(n1 = 0; n1 < num_beads-1; n1++) {
		for(n2 = n1+1; n2 < num_beads; n2++) 
    {
			q12mag = 0.0;
			for(j=0; j < DIMS; j++) {
				q12[j]=monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
				if((j==0 && wall_flag<3) || (j==2 && wall_flag<2) || (j==1 && wall_flag<1))
			    q12[j] = n_image(q12[j], maxsize[j]);
				q12mag += q12[j]*q12[j];
			}
			q12mag = sqrt(q12mag);

			if(sphere_pm->ev_type ==0) {
				for(i=0; i<DIMS; i++)
					force[i]=0.0;
			}
			/* WCA potential */
			else if(sphere_pm->ev_type ==1) {
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
			/* gaussian */
			else if(sphere_pm->ev_type ==2) {
				if(q12mag > sphere_pm->evcutoff*Rg)
					for(i=0; i<DIMS; i++)
						force[i]=0.0;
				else {
					aterm = aev/aev_alpha * exp(-q12mag*q12mag/aev_alpha);

					for(i=0; i<DIMS; i++) 
						force[i]= (-v*sphere_pm->kT*aterm*q12[i]);
				}	
			}

			for(i=0; i<DIMS; i++) {
				monomers[n1].force[i] -= force[i];
				monomers[n2].force[i] += force[i];
			}
      for(int m=0; m < DIMS; m++) {
        for(int n=0; n < DIMS; n++) {
          monomers[n1].stress_elas[m][n]-=monomers[n1].pos[m]*force[n];
          monomers[n2].stress_elas[m][n]+=monomers[n2].pos[m]*force[n];
        }
      }
		}
	}
}

/* Exc. Vol forces type = 0 HS, 1 LJ, 2 Gaussian */
//void ev_force_nlist(Float radius, struct sphere_param *sphere_pm, struct monomer *monomers)
//{
//	extern int max_x, max_y, max_z;
//  int maxsize[DIMS];
//  maxsize[0]=max_x;
//  maxsize[1]=max_y;
//  maxsize[2]=max_z;
//	int num_beads = sphere_pm->num_beads;
//  Float force[DIMS];
//	double Rg = sphere_pm->Ss;
//	double Nks =sphere_pm->nks;
//	double v = sphere_pm->sigma_k*sphere_pm->sigma_k*sphere_pm->sigma_k;
//	double aev_alpha = 4.0/3.0*Rg*Rg;
//	double aev = Nks*Nks*(1.0/((Pi*aev_alpha)*sqrt(Pi*aev_alpha)));
//	Float sigma = 2.0*radius;
//  sigma = 1/1.122;
//	Float cutoff = 1.122*sigma;
//	Float eps = sphere_pm->kT/sigma;
//  eps = sphere_pm->eps_WCA * sphere_pm->kT; // there should be a sigma in denominator?.
//
//  double sigma_attrac = sphere_pm->range_LJ;
//	double epsilon_LJ = sphere_pm->epsilon_LJ * sphere_pm->kT; 
//  double cutoff_LJ = sphere_pm->cutoff_LJ * sigma_attrac;
//
//  double r0 = sphere_pm->eq_morse;
//  double beta = sphere_pm->width_morse;
//  double epsMorse = sphere_pm->depth_morse;
//  double cutoff_morse = sphere_pm->cutoff_morse;
//
//  extern int wall_flag;
//  int i,j;
//  long n1, n2;
//  Float aterm;
//  Float q12[DIMS];
//  Float q12mag=0.0;
//  int sameParticle=0;
//  FILE *stream;
//
////	for(i=0; i < DIMS; i++) { force[i]=0.0; }
//	
//  for(n1=0; n1 < num_beads; n1++) {
//		for(n2=1; n2 <= monomers[n1].list[0]; n2++) 
//    {
//	    for(i=0; i < DIMS; i++) 
//        force[i]=0.0; 
//			// calculate the monomer-monomer distance
//      q12mag=0.0;
//			for(j=0; j<DIMS; j++) {
//				q12[j] = monomers[n1].pos_pbc[j] - monomers[monomers[n1].list[n2]].pos_pbc[j];
//				if((j==0 && wall_flag<3) || (j==2 && wall_flag<2) || (j==1 && wall_flag<1)) {
//					q12[j] = n_image(q12[j], maxsize[j]); 
//        }
//				q12mag += q12[j]*q12[j];
//			}	
//			q12mag = sqrt(q12mag);
//
//      // check whether monomers on the same particle
//      if(monomers[n1].sphere_id == monomers[monomers[n1].list[n2]].sphere_id)
//            sameParticle=1;
//      else  sameParticle = 0; 
//
//      if(sameParticle==1) 
//      {
//        if(sphere_pm->ev_type==0) 
//        {
//				  for(i=0; i<DIMS; i++) {
//					  force[j]=0.0; }
//			  }
//			  else if(sphere_pm->ev_type ==1)  // WCA potential
//        { 
//				  if(q12mag > cutoff) { 
//				    for(i=0; i<DIMS; i++) 
//						  force[i]=0.0; 
//          }
//				  else if(q12mag < cutoff) {
//					  aterm = (sigma/q12mag);
//					  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
//					  for(i=0; i<DIMS; i++) {
//						  force[i]=-24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag; }
//				  }
//			  }
//			  else if(sphere_pm->ev_type ==2) 
//        {
//				  aterm = aev/aev_alpha * exp(-q12mag*q12mag/aev_alpha);
//				  for(i=0; i<DIMS; i++) {
//					  force[i]= (-v*sphere_pm->kT*aterm*q12[i]); }
//			  }
//      }
//      else
//      {
//        if(sphere_pm->attrac_type == 1) 
//        {
//          if(q12mag > cutoff_LJ) {
//            for(i=0; i < DIMS; i++)
//              force[i] = 0.0;
//          }
//          else if(q12mag < cutoff_LJ) {
//            aterm = (sigma_attrac/q12mag);
//            aterm = aterm*aterm*aterm*aterm*aterm*aterm;
//            for(i=0; i<DIMS; i++) {
//              force[i]=-24.0*(epsilon_LJ*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag; }
//          }
//        }
//        if(sphere_pm->attrac_type == 2)
//        {
//          double dr = q12mag - r0;
//          double beta_dr = beta*dr;
//          if(q12mag >= cutoff_morse) {
//            for(i=0; i < DIMS; i++)
//              force[i] = 0.0;
//          }
//          else {
//            double mag = epsMorse*2*beta*(exp(-1.*beta_dr)  - exp(-2.*beta_dr));  // derivative of potential
//            for(i=0; i < DIMS; i++)
//              force[i] =  mag * q12[i] / q12mag;            
//          }
//        }      
//		  }
//			for(i=0; i<DIMS; i++) 
//      {
//			  monomers[n1].force[i] -= force[i];
//        monomers[n1].force_ev[i] -= force[i];
//				/*monomers[monomers[n1].list[n2]].force[i] += force[i]; */
//				monomers[n1].mmInteract[i] -= (force[i]*sphere_pm->monmass);  
//			}
//			//        for(i=0; i < DIMS; i++) { 
//			//          for(j=0; j < DIMS; j++) {
//			//            monomers[n1].stress[i][j] += ((-q12[i]) * monomers[n1].mmInteract[j]); // Zero it before using
//			//          } 
//			//        }
//		}
//	}
//}
void nonbonded_interaction(struct sphere_param *sphere_pm, struct monomer *monomers)
{
  int numBead = sphere_pm->num_beads;
  extern int wall_flag;
  extern int max_x, max_y, max_z;
  int maxsize[DIMS];
  maxsize[0] = max_x;
  maxsize[1] = max_y;
  maxsize[2] = max_z;
  double cutoff = 0.4;
  double sigma = cutoff / 1.122; 
  double eps = 2*sphere_pm->epsilon_WCA * sphere_pm->kT;

  for(int n1=0; n1 < numBead-1; n1++)  
  {
    for(int n2=n1+1; n2 < numBead; n2++)
    {
      int bonded=0;  
      for(int bond=1; bond <= monomers[n2].blist[0][0]; bond++) {
        if(n1 == monomers[n2].blist[bond][0])  
          bonded=1;
      }
      for(int bond=1; bond <= monomers[n1].blist[0][0]; bond++) {
        if(n2 == monomers[n1].blist[bond][0])  
          bonded=1;
      }

      double q12mag = 0;
      double q12[DIMS];
      double force[DIMS];
      double aterm = 0;
      if(bonded == 0)
      {
			  for(int j=0; j < DIMS; j++) {
				  q12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
				  if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
			      q12[j] = n_image(q12[j], maxsize[j]);
				  q12mag += q12[j]*q12[j];
			  }
			  q12mag = sqrt(q12mag);
				if(q12mag > cutoff) 
					for(int i=0; i < DIMS; i++)
						force[i]=0.0;      
				else {
          printf("nonbonded pair = (%d  %d)\n get_forces.c line 876\n", n1, n2);
					aterm = (sigma/q12mag);
					aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					for(int i=0; i < DIMS; i++) 
            // partial U partial r
						force[i] = -24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
				}
 			  for(int i=0; i < DIMS; i++) {
				  monomers[n1].force[i] -= force[i];
				  monomers[n2].force[i] += force[i];
			  }
      }
    }
  }
}

void ev_force_nlist(Float radius, struct sphere_param *sphere_pm, struct monomer *monomers)
{
	extern int max_x, max_y, max_z, wall_flag;
  int maxsize[DIMS];
  maxsize[0]=max_x;  maxsize[1]=max_y;  maxsize[2]=max_z;
  char fileName[100]; 
  FILE *stream;
  extern char *work_dir;

	double Rg = sphere_pm->Ss;
	double Nks = sphere_pm->nks;
	double v = sphere_pm->sigma_k*sphere_pm->sigma_k*sphere_pm->sigma_k;
	double aev_alpha = 4.0/3.0*Rg*Rg;
	double aev = Nks*Nks*(1.0/((Pi*aev_alpha)*sqrt(Pi*aev_alpha)));

	//Float sigma = 2.0*radius;
  //double cutoff = 1.122 * sigma;
  //Float eps = sphere_pm->kT/sigma; // there should be a sigma in denominator?
  //double sigma = 1 / 1.122;
  //double cutoff = 1.;
  extern double springLength_eq;
  springLength_eq=1.0;                    // Modification 20170315: temporary design

  // WCA potential
  double eps = sphere_pm->epsilon_WCA * sphere_pm->kT;
  //double sigma = springLength_eq / pow(2.,1./6.);
  double sigma = sphere_pm->eq_WCA / pow(2.,1./6.);
  //double cutoff = springLength_eq;
  double cutoff = sphere_pm->eq_WCA;
  // L-J potential
	double epsilon_LJ = sphere_pm->epsilon_LJ * sphere_pm->kT; 
  double sigma_attrac = sphere_pm->eq_LJ / pow(2., 1./6.);
  double cutoff_LJ = sphere_pm->cutoff_LJ * sigma_attrac;
  // Morse potential
  double r0 = sphere_pm->eq_morse;
  double beta = sphere_pm->width_morse;
  double epsMorse = sphere_pm->depth_morse * sphere_pm->kT;
  double cutoff_morse = sphere_pm->cutoff_morse;
  // WCA potential for initilization procedure
  double eps_init=0.005652;
  double sigma_init=1/pow(2.,1./6.);
  double cutoff_init=1.;

  double forceA=0., forceR=0., forceR_morse=0., forceR_wca=0.;
  int counterA=0, counterR=0, counterR_morse=0, counterR_wca=0;  
	int num_beads = sphere_pm->num_beads;
  double shortestDis[num_beads];
 
  for(int n1=0; n1 < num_beads; n1++) 
  {
    shortestDis[n1]=100.; 
		for(int neighbor=1; neighbor <= monomers[n1].list[0]; neighbor++) 
    {
      int n2 = monomers[n1].list[neighbor];
      if(monomers[n1].sphere_id != monomers[n2].sphere_id)
      {
        double q12[3], q12mag=0., force[3]={0.};    
			  // calculate the monomer-monomer distance
			  for(int j=0; j < DIMS; j++) {
				  q12[j] = monomers[n1].pos_pbc[j] - monomers[monomers[n1].list[neighbor]].pos_pbc[j];
				  if((j==0 && wall_flag<3) || (j==2 && wall_flag<2) || (j==1 && wall_flag<1)) {
					  q12[j] = n_image(q12[j], maxsize[j]); 
          }
				  q12mag += q12[j]*q12[j];
			  }	
			  q12mag = sqrt(q12mag);
        if(q12mag<shortestDis[n1])  shortestDis[n1]=q12mag;

        if(sphere_pm->attrac_type == 0) // WCA; For generating a initial configuration
        { 
				  //if(q12mag >= cutoff_init) { 
					//	  force[0]=0.;  force[1]=0.;  force[2]=0.;   
          //}
				  //else //else if(q12mag < cutoff_init) 
          //{
					//  double aterm = sigma_init / q12mag;
					//  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					//  for(int i=0; i < DIMS; i++) 
					//	  force[i] = -24.0*(eps_init*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
				  //}
				  if(q12mag < cutoff_init) 
          {
					  double aterm = sigma_init / q12mag;
					  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					  for(int i=0; i < DIMS; i++) 
						  force[i] = 24.0*(eps_init*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
				  }
				  //else { 
					//	  force[0]=0.;  force[1]=0.;  force[2]=0.;   
          //}
			  }
        else if(sphere_pm->attrac_type == 1) // LJ potential
        {
          if(q12mag < cutoff_LJ) 
          {
            double aterm = sigma_attrac / q12mag;
            aterm = aterm*aterm*aterm*aterm*aterm*aterm;
            for(int i=0; i<DIMS; i++) {
              force[i] = 24.0*(epsilon_LJ*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
              monomers[n1].force_interA[i] = force[i];            
            }
          }
        }
        else if(sphere_pm->attrac_type == 2) // Morse potential
        {
          double dr = q12mag - r0;
          double beta_dr = beta*dr;
          //if(q12mag >= cutoff_morse) {
          //  force[0] = 0.;  force[1]=0.; force[2]=0.;
          //}
          //else {
          //  double mag = epsMorse*2*beta*(exp(-1.*beta_dr) - exp(-2.*beta_dr)); // derivative of potential
          //  force[0] =  mag * q12[0] / q12mag;   
          //  force[1] =  mag * q12[1] / q12mag; 
          //  force[2] =  mag * q12[2] / q12mag;     
          //}
          if(q12mag < cutoff_morse)
          {
            double mag = epsMorse*2*beta*(exp(-1.*beta_dr) - exp(-2.*beta_dr)); // derivative of potential
            force[0] =  -mag * q12[0] / q12mag;   
            force[1] =  -mag * q12[1] / q12mag; 
            force[2] =  -mag * q12[2] / q12mag;
            monomers[n1].force_interA[0]=force[0];
            monomers[n1].force_interA[1]=force[1];
            monomers[n1].force_interA[2]=force[2];    
          }
          //else {
          //  force[0]=0.;  force[1]=0.;  force[2]=0.;
          //}
        }
        else if(sphere_pm->attrac_type ==3) // blend of LJ and Morse potentials
        {                                   // Note 20170814: Revision needed!
          double dr = q12mag - r0;
          double beta_dr = beta*dr;
          if(q12mag > cutoff_morse) {
            for(int i=0; i < DIMS; i++)
              force[i] = 0.;
          }
          else if(q12mag <= cutoff_morse && q12mag > r0) {
            double mag = epsMorse*2*beta*(exp(-1.*beta_dr) - exp(-2.*beta_dr));  // derivative of potential
            for(int i=0; i < DIMS; i++)
              force[i] =  mag * q12[i] / q12mag;
          }
          else if(q12mag <=  sphere_pm->eq_LJ/*sigma_attrac*pow(2., 1./6.)*/) {
					  double aterm = (sigma_attrac/q12mag);
					  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					  for(int i=0; i < DIMS; i++) 
						  force[i]=-24.0*(epsilon_LJ*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
          }
        }
        else if(sphere_pm->attrac_type == 4)  // blend of WCA and Morse potentials
        {
          double dr = q12mag - r0;
          double beta_dr = beta*dr;
          //if(q12mag >= cutoff_morse) {
          //  force[0]+=0.;  force[1]+=0.; force[2]+=0.;
          //}
          //else {
          //  double mag = epsMorse*2*beta*(exp(-1.*beta_dr) - exp(-2.*beta_dr)); // derivative of potential
          //  force[0] +=  mag * q12[0] / q12mag;   
          //  force[1] +=  mag * q12[1] / q12mag; 
          //  force[2] +=  mag * q12[2] / q12mag;     
          //}
				  //if(q12mag >= cutoff) { 
					//	  force[0]+=0.;  force[1]+=0.;  force[2]+=0.;   
          //}
				  //else //else if(q12mag < cutoff) 
          //{
					//  double aterm = sigma/q12mag;
					//  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					//  for(int i=0; i < DIMS; i++) 
					//	  force[i] -= 24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
				  //}
          if(q12mag < cutoff_morse) {
            double mag = epsMorse*2*beta*(exp(-1.*beta_dr) - exp(-2.*beta_dr)); // derivative of potential
            double forceTemp[3];
            forceTemp[0] = -mag * q12[0] / q12mag;
            forceTemp[1] = -mag * q12[1] / q12mag;
            forceTemp[2] = -mag * q12[2] / q12mag;
            force[0] +=  forceTemp[0];   
            force[1] +=  forceTemp[1]; 
            force[2] +=  forceTemp[2];
            monomers[n1].force_interA[0]=forceTemp[0];
            monomers[n1].force_interA[1]=forceTemp[1];
            monomers[n1].force_interA[2]=forceTemp[2];

            // Record force magnitude
            if(q12mag < r0) {
              double magnitude = sqrt(iproduct(forceTemp,forceTemp));
              forceR_morse += magnitude;
              counterR_morse++;
              forceR += magnitude; 
            }
            else {
              forceA += sqrt(iproduct(forceTemp,forceTemp));
              counterA++; 
            }
          }
          //else {
          //  force[0]+=0.;  force[1]+=0.; force[2]+=0.;
          //}
				  if(q12mag < cutoff) 
          {
					  double aterm = sigma/q12mag;
					  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
            double forceTemp[3];
            forceTemp[0]=24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[0]/q12mag/q12mag;
            forceTemp[1]=24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[1]/q12mag/q12mag;
            forceTemp[2]=24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[2]/q12mag/q12mag;
						force[0] += forceTemp[0];
            force[1] += forceTemp[1];
            force[2] += forceTemp[2];  
            monomers[n1].force_interR[0]=forceTemp[0];
            monomers[n1].force_interR[1]=forceTemp[1];
            monomers[n1].force_interR[2]=forceTemp[2];

            // Record force magnitude
            double magnitude=sqrt(iproduct(forceTemp,forceTemp));
            forceR_wca += magnitude;
            counterR_wca++;
            forceR += magnitude;
            counterR++;
				  }
				  //else { 
					//	  force[0]+=0.;  force[1]+=0.;  force[2]+=0.;   
          //}
        }
        else if(sphere_pm->attrac_type == 5)  // repulsive WCA potential
        { 
				  //if(q12mag >= cutoff) { 
					//	  force[0]=0.;  force[1]=0.;  force[2]=0.;   
          //}
				  //else //else if(q12mag < cutoff) 
          //{
					//  double aterm = sigma / q12mag;
					//  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					//  for(int i=0; i < DIMS; i++) 
					//	  force[i] = -24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
				  //}
				  if(q12mag < cutoff) 
          {
					  double aterm = sigma / q12mag;
					  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					  for(int i=0; i < DIMS; i++) {
						  force[i] = 24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
              monomers[n1].force_interR[i]=force[i];   
            }
            double magnitude = sqrt(iproduct(force,force));
            forceR += magnitude;
            counterR++;
				  }
				  //else { 
					//	  force[0]=0.;  force[1]=0.;  force[2]=0.;   
          //}
        }
        monomers[n1].force[0]+=force[0];  
        monomers[n1].force[1]+=force[1]; 
        monomers[n1].force[2]+=force[2]; 
        monomers[n1].force_inter[0]=force[0];  
        monomers[n1].force_inter[1]=force[1];  
        monomers[n1].force_inter[2]=force[2];
        for(int m=0; m < DIMS; m++) { // Calculate particle stress
          for(int n=0; n < DIMS; n++) { 
            monomers[n1].stress_int_v1[m][n] += 0.5*q12[m]*force[n];
            monomers[n1].stress_int_v2[m][n] -= 0.5*q12[m]*force[n];
          }
        } 
        //double forcemag= sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]); 
        //if(q12mag < 0.4) 
        //{
        //  sprintf(fileName,"%s/data/warningMes.dat",work_dir);
        //  stream=fopen(fileName,"a");
        //  fprintf(stream,"vertex-vertex distance = %f  force = %f\n", q12mag,forcemag);
        //  fclose(stream);
		    //}     
      }
	  }
  }

  // Record force magnitudes
  if(forceR>0)  forceR /= counterR;
  if(forceR_wca>0)  forceR_wca /= counterR_wca;
  if(forceR_morse>0) forceR_morse /= counterR_morse;
  if(forceA>0) forceA /= counterA; 
  if(counterR>0 || counterA>0) {  // Modification 20170811
    sprintf(fileName,"%s/data/interPartForce.dat",work_dir);
    stream=fopen(fileName,"a");
    fprintf(stream,"%f    %f    %f    %f\n",forceA,forceR,forceR_wca,forceR_morse);
    fclose(stream);
  }

  // write shortest distance between vertices on different particles 
  sprintf(fileName,"%s/data/shortestDis.dat",work_dir);
  stream=fopen(fileName,"a");
  for(int n=0; n<num_beads; n++)
    fprintf(stream,"%f ",shortestDis[n]);
  fprintf(stream,"\n");
  fclose(stream);
}

void nonbonded_interaction_nlist(struct sphere_param *sphere_pm, struct monomer *monomers)
{
  int num_beads = sphere_pm->num_beads;
	extern int max_x, max_y, max_z, wall_flag;
  int maxsize[DIMS];
  maxsize[0] = max_x;
  maxsize[1] = max_y;
  maxsize[2] = max_z;
  double cutoff = 0.7;
  double sigma = cutoff / 1.122; 
	double eps = 0.005652; //sphere_pm->epsilon_WCA * sphere_pm->kT;  // Modification 20170802: Temporary!!!
	FILE *stream;
  char fileName[100];
  extern char *work_dir;

  for(int n1=0; n1 < num_beads; n1++) {
		for(int listLabel=1; listLabel <= monomers[n1].list[0]; listLabel++) 
    {
      int n2 = monomers[n1].list[listLabel];
      int bonded = 0;
      if(monomers[n1].sphere_id == monomers[n2].sphere_id) 
      { 
        for(int bond=1; bond <= monomers[n2].blist[0][0]; bond++) {
          if(n1 == monomers[n2].blist[bond][0])  
            bonded=1;
        }
        for(int bond=1; bond <= monomers[n1].blist[0][0]; bond++) {
          if(n2 == monomers[n1].blist[bond][0])  
            bonded=1;
        }
        double q12mag = 0;
        double q12[DIMS];
        double force[DIMS]={0.};
        double aterm = 0;
        if(bonded == 0)
        {
			    for(int j=0; j < DIMS; j++) {
				    q12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
				    if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
			        q12[j] = n_image(q12[j], maxsize[j]);
				    q12mag += q12[j]*q12[j];
			    }  
			    q12mag = sqrt(q12mag);
				  if(q12mag > cutoff) {
						force[0]=0.;  force[1]=0.;  force[2]=0.;      
          }
				  else {              
					  aterm = (sigma/q12mag);
					  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
					  for(int i=0; i < DIMS; i++) 
              // partial U partial r
						  force[i] = -24.0*(eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;

            //sprintf(fileName,"%s/data/warningMes.dat",work_dir);
            //stream=fopen(fileName,"a");
            //fprintf(stream,"dimple distance = %f  force = %f\n", q12mag, 
            //        sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]));
            //fclose(stream);

				  }
 			    
				  monomers[n1].force[0] -= force[0]; monomers[n1].force[1] -= force[1]; 
          monomers[n1].force[2] -= force[2];
         
          monomers[n1].force_nonbonded_intra[0] -= force[0];
	        monomers[n1].force_nonbonded_intra[1] -= force[1];	
          monomers[n1].force_nonbonded_intra[2] -= force[2];
          //if(q12mag < 0.2)
          //{
            //sprintf(fileName,"%s/data/warningMes.dat",work_dir);
            //stream=fopen(fileName,"a");
            //fprintf(stream,"dimple distance = %f  force = %f\n", q12mag, 
            //        sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]));
            //fclose(stream);           
          //}
        }
      }
    }
  }
}

/* Calculate spring force : type = 0 FENE, 1 WLC */
void wall_force(struct sphere_param *sphere_pm, struct monomer *monomers)
{
	extern int n_proc, num_x, max_x, max_y, max_z, wall_flag;
	long n1, n2;
	int i,j;
	int num_beads = sphere_pm->num_beads;
	Float r1;
	Float mag;
	Float wally, wallz;
	Float wall_aev;
	Float force[DIMS];
	double maxsize[DIMS], monpos[DIMS];
	FILE *stream;
  char fileName[100];
  extern char* work_dir;
	maxsize[0]=(double)max_x;
	maxsize[1]=(double)max_y;
	maxsize[2]=(double)max_z;

  double forceWall=0.;  // Modification 20170811
  int counterWall=0;    // Modification 20170811

	if(monomers[0].radius*3 > 1.0) 
  {
		//    wally = (maxsize[1]-1.0)/2.0-monomers[0].radius;
		//    wallz = (maxsize[2]-1.0)/2.0-monomers[0].radius;

    /* Note 20170314:
    -----------------    
    The origin is at the center of the box. 
    In the particle coordinate system, the endpoints of the fluid domain are 0 and max-1. 
    Shifting the origin to the center of box, the endpoints of the fluid domain in the 
    particle coordinate system are (max-1)/2 and -(max-1)/2. 
    1.2 is an arbitrary factor! 
    */
		wally = (maxsize[1]-1.0)/2.0 - 1.2;
		wallz = (maxsize[2]-1.0)/2.0 - 1.2;

    /* Modification 20170710:
    ------------------------- 
    Temporary modification! Now kT=1.
    */   
		//wall_aev = 25.0*sphere_pm->kT / (sphere_pm->sigma_k*sphere_pm->sigma_k*sphere_pm->sigma_k);    
    wall_aev = 14.8;    
	}
	else 
  {
		wally = (maxsize[1]-1.0)/2.0-1.0;
		wallz = (maxsize[2]-1.0)/2.0-1.0;
		wall_aev = 25.0*sphere_pm->kT/(sphere_pm->sigma_k*sphere_pm->sigma_k*sphere_pm->sigma_k);
	}

	for(n1=0; n1 < num_beads; n1++) 
  {
		if((int)(monomers[n1].pos[0]/num_x) == n_proc)
			monpos[0] = monomers[n1].pos[0] - ((maxsize[0]-1.0)/2.0+1.0);
		for(j=1; j<DIMS; j++)
			monpos[j] = monomers[n1].pos[j] - ((maxsize[j]-1.0)/2.0);

		if(wall_flag >= 1) 
    {
			if(monpos[1] < -wally)
				force[1] = wall_aev * (monpos[1]+wally)*(monpos[1]+wally);
			else if(monpos[1] > wally)
				force[1] = -wall_aev * (monpos[1]-wally)*(monpos[1]-wally);
			else
				force[1]=0.0;

			monomers[n1].force[1] += force[1];
      monomers[n1].force_wall[1] += force[1];
      for(int m=0; m < DIMS; m++) 
        monomers[n1].stress_wall[m][1] += monomers[n1].pos[m]*force[1];
      
      forceWall += force[1];  // Modification 20170811
      counterWall++;          // Modification 20170811
 
			if(fabs(monpos[1]) > (maxsize[1]-1.0)/2.0) 
      {
        sprintf(fileName,"%s/data/warningMes.dat",work_dir);
        stream=fopen(fileName,"a");
				fprintf(stream,"bead %d out of ybounds, pos=%f %f %f \n",
        n1, monomers[n1].pos[0], monomers[n1].pos[1], monomers[n1].pos[2]);
			  fprintf(stream,"fluid velocity = %f %f %f\nvel=%f %f %f\n",
        monomers[n1].fluid_vel[0], monomers[n1].fluid_vel[1],monomers[n1].fluid_vel[2],
        monomers[n1].vel[0], monomers[n1].vel[1], monomers[n1].vel[2]);
				fprintf(stream,"force = %f %f %f\ninternal force=%f %f %f\n", 
        monomers[n1].force[0], monomers[n1].force[1], monomers[n1].force[2], 
        monomers[n1].f_int[0], monomers[n1].f_int[1], monomers[n1].f_int[2]);
				fprintf(stream,"wall force = %f\n", force[1]);
        fclose(stream);
				exit(18);
			}
		}
		if(wall_flag == 2) 
    {
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

      monomers[n1].force_wall[1] += force[1];
      monomers[n1].force_wall[2] += force[2];

			if(fabs(monpos[1]) > (maxsize[1]-1.0)/2.0 || fabs(monpos[2]) > (maxsize[2]-1.0)/2.0) 
      {
				printf(" bead %d out of yzbounds, pos=%le %le %le \n", n1, 
               monomers[n1].pos[0], monomers[n1].pos[1], monomers[n1].pos[2]);
				printf(" fluid velocity = %le %le %le  vel=%le %le %le\n", monomers[n1].fluid_vel[0], 
               monomers[n1].fluid_vel[1], monomers[n1].fluid_vel[2], monomers[n1].vel[0], 
               monomers[n1].vel[1], monomers[n1].vel[2]);
				printf(" force = %le %le %le internal force=%le %le %le\n", monomers[n1].force[0], 
               monomers[n1].force[1], monomers[n1].force[2], monomers[n1].f_int[0], 
               monomers[n1].f_int[1], monomers[n1].f_int[2]);
				printf(" wall force = %le %le\n", force[1], force[2]);
				exit(18);
			}
		}
	}

  forceWall /= counterWall;  // Modification 20170811
  if(forceWall>0) {        // Modification 20170811
    sprintf(fileName,"%s/data/wallForce.dat",work_dir);
    stream=fopen(fileName,"a");
    fprintf(stream,"%f\n",forceWall); 
    fclose(stream);
  }
}

/* Hydrodynamic forces */
void hi_force(struct monomer *mon, Float ***velcs_df, double fric, double tau, double dt, 
     int num_beads, Float mon_mass, VSLStreamStatePtr rngstream)
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
	int x[DIMS], x1[DIMS];
	int qxp, qyp, qzp;
	int qxm, qym, qzm;
	int maxsize[DIMS];
	double x2[DIMS];
	double **vel, *tmp_pp;
	double vel_node[8][DIMS], rho_node[8];
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
	for(n1=0; n1 < num_beads; n1++) 
  {
		for(d=0; d<DIMS; d++)
			x[d] = mon[n1].pos[d];  // x is an integer
		fluid_vel(n1, vel, x, velcs_df, mon);
		//    if(MD_tau == tau)
		for(d=0; d<DIMS; d++) 
			mon[n1].fluid_vel[d]=vel[n1][d];

		for(d=0; d<DIMS; d++) {
			mon[n1].fricforce[d] = -fric * (mon[n1].vel[d]-vel[n1][d]);
			//mon[n1].fricforce[d]=-fric*(mon[n1].vel[d]);
			mon[n1].drag[d] = -fric * Rho_Fl * (mon[n1].vel[d] - vel[n1][d]); 
      //mon[n1].v_diff[d] = mon[n1].vel[d] - vel[n1][d]; //20170125
		}
	}
	/* Generate the random force */
	if(add_noise >= 1 && add_noise !=3)
		fluc_force(fforce, num_beads*DIMS, rngstream);

	/* redistribute momentum to the surrounding lattice */
	n1=0;
	while(n1<num_beads) 
  {
		for(d=0; d<DIMS; d++) {
			x[d] = (int)floor(mon[n1].pos[d]); /* 0 <= floor(monomer position) <= max_size-1 */
			x2[d]=mon[n1].pos[d]-(double)x[d]; /* translate the coordinates of the monomer */
		}
		mon_proc = box(x[0],maxsize[0]) / num_x;
		if(n_proc == mon_proc) 
    {
			for(i=0; i<DIMS; i++) {
				force[i] = Rho_Fl*mon[n1].fricforce[i];
				force[i] += fforce[n1*DIMS+i];
			}
			for(i=0; i<DIMS; i++) {
				mon[n1].f_fluc[i]=fforce[n1*DIMS+i];
        mon[n1].force_fluc[i] = fforce[n1*DIMS+i];
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

			for(k=0; k<=1; k++) 
      {
				x1[2] = box(x[2]+k, maxsize[2])+1.0;
				if(k==0)    temp[2] = (1-x2[2]);
				else  			temp[2] = x2[2];
				for(j=0; j<=1; j++) 
        {
					x1[1] = box(x[1]+j, maxsize[1])+1.0;
					if(j==0)    temp[1] = (1-x2[1]);
					else  			temp[1] = x2[1];
					for(i=0; i<=1; i++) 
          { 
						x1[0] = box((x[0]+i), maxsize[0])%num_x+1.0;
						if(i==0)    temp[0] = 1-x2[0];
						else  			temp[0] = x2[0];
						nxy = x1[0]*(max_y+2) + x1[1];
						for(q=0; q < Num_Dir; q++) 
							velcs_df[nxy][q][x1[2]] += dmomentum[q]*temp[0]*temp[1]*temp[2]; 

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

void fluid_vel(int n1, double **vel, int x[DIMS], Float ***velcs_df, struct monomer *mon)
{
	extern int wall_flag;
	extern int max_x, max_y, max_z;
	extern int num_x;
	int i,j,k;
	int d,q,n,nxy;
	int x1[DIMS], maxsize[DIMS];
	double x2[DIMS];
	double vel_node[8][DIMS], rho_node[8];
	double temp[DIMS];
	double rho;
	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;

	for(i=0; i < 8; i++)
		rho_node[i]=0.0;
  // x1 is an integer served as the Eulerian coordinate
	for(k=0; k<=1; k++) {  
		x1[2] = box(x[2]+k, maxsize[2]) + 1;
		for(j=0; j<=1; j++) {
			x1[1] = box(x[1]+j, maxsize[1]) + 1;
			for(i=0; i<=1; i++) {
				x1[0] = box((x[0]+i), maxsize[0])%num_x + 1.0;
				nxy = x1[0]*(maxsize[1]+2) + (x1[1]);
				n = i*4+j*2+k;
				rho_node[n] = 0.0;
			  vel_node[n][0] = 0.0;
        vel_node[n][1] = 0.0;
        vel_node[n][2] = 0.0;
				for(q=0; q < Num_Dir; q++) { 
					rho_node[n]   += velcs_df[nxy][q][ x1[2] ];
					vel_node[n][0]+= c_x[q]*velcs_df[nxy][q][x1[2]];
					vel_node[n][1]+= c_y[q]*velcs_df[nxy][q][x1[2]];
					vel_node[n][2]+= c_z[q]*velcs_df[nxy][q][x1[2]];
				}
        //printf("pre node rho %le (%d %d %d) vel (%le %le %le)\n", rho_node[n], x1[0], 
        //       x1[1], x1[2], vel_node[n][0], vel_node[n][1], vel_node[n][2]);
			}
		}
	}
//modified 20160919
	/* wall boundary conditions */
//	if(wall_flag >= 1) { 
//		if(x[1] == 0) 
//			for(i=0; i<=1; i++) 
//				for(k=0; k<=1; k++) 
//					for(d=0; d<DIMS; d++)
//						vel_node[i*4+k][d]=-vel_node[i*4+2+k][d];
//
//		if(x[1] == maxsize[1]-1) 
//			for(i=0; i<=1; i++)
//				for(k=0; k<=1; k++) 
//					for(d=0; d<DIMS; d++)
//						vel_node[i*4+2+k][d]=-vel_node[i*4+k][d];
//	}
//
//	if(wall_flag == 2) { 
//		if(x[2] == 0) 
//			for(i=0; i<=1; i++) 
//				for(j=0; j<=1; j++) 
//					for(d=0; d<DIMS; d++)
//						vel_node[i*4+j*2][d]=-vel_node[i*4+j*2+1][d];
//
//		if(x[2] == maxsize[2]-1) 
//			for(i=0; i<=1; i++)
//				for(j=0; j<=1; j++) 
//					for(d=0; d<DIMS; d++)
//						vel_node[i*4+j*2+1][d]=-vel_node[i*4+j*2][d];
//	}

	/* Use trilinear interpolation to get the fluid velocity at the monomer position */
	for(d=0; d<DIMS; d++) {
		vel[n1][d] = 0.0;
		x2[d] = mon[n1].pos[d]-(double)x[d];   /* translate the coordinates of the monomer */
	}
	for(k=0; k<=1; k++) {
		if(k==0)    temp[2] = 1.0-x2[2];
		else        temp[2] = x2[2];
		for(j=0; j<=1; j++) {
			if(j==0)    temp[1]=1.0-x2[1];
			else 				temp[1]=x2[1];
			for(i=0; i<=1; i++) { 
				if(i==0)    temp[0] = 1.0-x2[0];
				else        temp[0] = x2[0];
				n = 4*i+2*j+k;
				for(d=0; d < DIMS; d++) 
					vel[n1][d] += vel_node[n][d]*temp[0]*temp[1]*temp[2];
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

void ParticleStress(struct sphere_param *sphere_pm, struct monomer *monomers)
{
  /* This function only calculate particle stress attributing to particle 
  elasticity, bending rigidity, volume and area constraints, and wall force. 
  The extra stress caused by inter-particle interaction is calculated in the 
  ev_force function. To obtain results consistent with Kruger's, the 
  frictional force is excluded when particle stress is evaluated. */
	
	int numBead = sphere_pm->num_beads;
	for(int n=0; n < numBead; n++) {
		for(int i=0; i < DIMS; i++) {
			for(int j=0; j < DIMS; j++) {
				monomers[n].stress[i][j] += -1. * monomers[n].pos[i] *
        (monomers[n].force[j]- monomers[n].force_pp[j] - monomers[n].drag[j]);
			}
		}
	}
}

void spreading(struct vector **forceDen, struct monomer *mon, int num_beads) 
{
  extern int max_x, max_y, max_z;
  int    x[DIMS];
  double x2[DIMS];
  int    x1[DIMS];
  int maxsize[DIMS];
  double temp[DIMS];
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  // Zero the external force density
  for(int x=1; x <= max_x; x++) {
    for(int y=1; y <= max_y; y++) {
      int xy = x*(max_y+2) + y;
      for(int z=1; z <= max_z; z++) {
        forceDen[xy][z].x = 0;
        forceDen[xy][z].y = 0;
        forceDen[xy][z].z = 0;
      }
    }
  }
  
  for(int n=0; n < num_beads; n++)
  {
    for(int d=0; d < DIMS; d++) {
      x[d] = (int)floor(mon[n].pos[d]);
      x2[d] = mon[n].pos[d]-(double)x[d];  // The fractional part of x[d]
    }   
    // weight factors
    for(int k=0; k <= 1; k++)
    {
      x1[2] = box(x[2]+k, maxsize[2]) + 1.0;
      if(k == 0)    temp[2] = 1-x2[2]; 
      else          temp[2] = x2[2];     
      for(int j=0; j <= 1; j++)
      {
        x1[1] = box(x[1]+j, maxsize[1]) + 1.0;
        if(j == 0)    temp[1] = 1-x2[1]; 
        else          temp[1] = x2[1];     
        for(int i=0; i <= 1; i++)
        {
          x1[0] = box(x[0]+i, maxsize[0]) + 1.0;
          if(i == 0)    temp[0] = 1-x2[0]; 
          else          temp[0] = x2[0];   
          int nxy = x1[0]*(max_y+2) + x1[1];
          forceDen[nxy][x1[2]].x += mon[n].force[0]*temp[0]*temp[1]*temp[2];
          forceDen[nxy][x1[2]].y += mon[n].force[1]*temp[0]*temp[1]*temp[2];
          forceDen[nxy][x1[2]].z += mon[n].force[2]*temp[0]*temp[1]*temp[2];
//printf("forceDen[%d][%d]=(%f, %f, %f)\n",nxy,x1[2],forceDen[nxy][x1[2]].x,forceDen[nxy][x1[2]].y,forceDen[nxy][x1[2]].z);
        }
      }
    }
  }
//PAUSE
}

void interpolation(int num_beads, Float ***velcs_df, struct vector **forceDen, double dt, 
     struct monomer *mon)
{
  extern int max_x, max_y, max_z;
  int lagrange[DIMS], euler[DIMS];
  double x2[DIMS];
  double vel_node[8][DIMS], rho_node[8];
  double temp[DIMS];
  int maxsize[DIMS];
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  for(int n=0; n < num_beads; n++)  // loop over all beads
  {
    // Find the integral part of the lagrangian coordinate of a monomer
    lagrange[0] = (int)floor(mon[n].pos[0]);
    lagrange[1] = (int)floor(mon[n].pos[1]);
    lagrange[2] = (int)floor(mon[n].pos[2]);

    // Find the fluid nodes surrounding a monomer and calculate their density and velocity
    for(int k=0; k <= 1; k++) {
      euler[2] = box(lagrange[2] + k, maxsize[2]) + 1;
      for(int j=0; j <= 1; j++) {
        euler[1] = box(lagrange[1] + j, maxsize[1]) + 1;
        for(int i=0; i <= 1; i++) {
          euler[0] = box(lagrange[0] + i, maxsize[0]) + 1;
          int nxy = euler[0] * (maxsize[1]+2) + (euler[1]);
          int index = i*4 + j*2 +k;
          rho_node[index] = 0.0;     // Zero the density of local fluid surrounding the monomer
          vel_node[index][0] = 0.0;  // Zero the velocity of local fluid surrounding the monomer    
          vel_node[index][1] = 0.0;
          vel_node[index][2] = 0.0;    
          for(int q=0; q < Num_Dir; q++) {
            rho_node[index]    += velcs_df[nxy][q][euler[2]];
            vel_node[index][0] += c_x[q]*velcs_df[nxy][q][euler[2]];
            vel_node[index][1] += c_y[q]*velcs_df[nxy][q][euler[2]];
            vel_node[index][2] += c_z[q]*velcs_df[nxy][q][euler[2]];
          }
//          rho_node[index] += Rho_Fl;  
          vel_node[index][0] /= rho_node[index];
          vel_node[index][1] /= rho_node[index];
          vel_node[index][2] /= rho_node[index];

          vel_node[index][0] += (0.5*forceDen[nxy][euler[2]].x*dt/rho_node[index]);
          vel_node[index][1] += (0.5*forceDen[nxy][euler[2]].y*dt/rho_node[index]);
          vel_node[index][2] += (0.5*forceDen[nxy][euler[2]].z*dt/rho_node[index]); 
        }
      }
    }
    // Use a two-point linear interpolation function to evaluate the velocity of a monomer
    for(int d=0; d < DIMS; d++) {
      mon[n].vel[d] = 0.0;  // Should I zero it here??
      x2[d] = mon[n].pos[d] - (double) lagrange[d];  // The fractional part
    }
    for(int k=0; k <= 1; k++)
    {
      if(k==0) temp[2] = 1.0 - x2[2];
      else     temp[2] = x2[2];
      for(int j=0; j <= 1; j++)
      {
        if(j==0)  temp[1] = 1.0 - x2[1];
        else    temp[1] = x2[1];
        for(int i=0; i <= 1; i++)
        {
          if(i==0)  temp[0] = 1.0 - x2[0];
          else       temp[0] = x2[0];
          int index = i*4 + j*2 + k;
          for(int d=0; d < DIMS; d++)
            mon[n].vel[d] += (vel_node[index][d]*temp[0]*temp[1]*temp[2]);
        }
      }
    }
  }
}


