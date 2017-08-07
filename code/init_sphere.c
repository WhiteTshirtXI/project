/* Model sphere as an icosahedron w/ vertices joined by springs */
/* Model RBC as compressed level 2 spheres from template */      
/* spheres always inserted before RBC */                                            
/* Coupled with LBE code to add hydrodynamic forces */ 

#include "header.h"
#include <time.h>
#define NVERTICES 12   /* number of vertices for an icosahedron */                    
#define NBONDS 5       /* number of bonds for a vertex on an icosahedron */          
#define MRBC 162       /* number of vertices for RBC */                              

void setup(double ***v, int blist[NVERTICES][NBONDS+1][3], double radius[NTYPES]);
void setup_template(double **, int ***, int, double, char *);
int check_bond(int k, int m, struct monomer *monomers);
int check_face(int, int, int, struct face *, int);
int get_midpoint(int a, int b, struct monomer *monomers, int **midpoints);

int sphere_init(struct sphere_param *sphere_pm, struct sphere *spheres, struct monomer 
        *monomers, struct face *faces, char *work_dir)
{
	extern int n_proc;
	extern int max_x, max_y, max_z;
	extern int wall_flag;
	int i,j,k,m,f, pnum, m_temp;
	int n, level;
	int n1, n2, n3, n4;
	int n12, n13, n14, n23, n24;
	int nx, ny, nz;
	int trial_count=0;
	int rot_flag;
	int temp;
	int type, type1, start, start0;
	int Nsphere = sphere_pm->Nsphere;
	int num_beads = sphere_pm->num_beads;
	int maxNbeads;
	int Ntype[NTYPES];
	int nlevel[NTYPES];
	int Nbeads[NTYPES];
	int Nfaces[NTYPES];
	int blist_ico[NVERTICES][NBONDS+1][3]; /* blist for the starting icosahedron */
	int bcount;
	int ***blist_temp, ***blist_rbc;
	int **midpoints;
	double maxsize[DIMS];
	double bondlength, radius[NTYPES];
	double Vfactor[NTYPES];
	double var;
	double dx, dy, dz, r2;
	double **centers;                              /* center position of spheres */
	double theta, phi;
	double axis[DIMS][DIMS];                             /* x,y,z axis of sphere */
	double ***v;                                                /* the bond list */
	double p_mid[DIMS], dr[DIMS];
	char filename[200];
	FILE *stream=0;
	int overlap = FALSE;
	int test=TRUE;                              // It's always TRUE.
	char test_file[200] = {"data/init.dat"};
	FILE *file_ptr=0;
  int trial_count_max=150000;
  
  int overlap_pp=0;  // Modification 20170320: for recording the frequency of overlap.
  int overlap_pw=0;
  int rotation=0;

	srand((unsigned)time(0));
	maxsize[0] = (double)max_x;
	maxsize[1] = (double)max_y;
	maxsize[2] = (double)max_z;

	for (i=0 ; i < NTYPES ; i++) 
	{
		Ntype[i]  = sphere_pm->Ntype[i];
		nlevel[i] = sphere_pm->nlevel[i];
		Nbeads[i] = sphere_pm->N_per_sphere[i];
		Nfaces[i] = sphere_pm->face_per_sphere[i];
	}
	// Memory allocation 
	centers = (double**) calloc(Nsphere, sizeof(double*));
	if (centers == 0)  fatal_err ("cannot allocate centers", -1);

	midpoints = (int **) calloc(num_beads, sizeof(int *));
	if (midpoints == 0) fatal_err("cannot allocate midpoints", -1);

	v = (double ***) calloc(NTYPES, sizeof(double **));
	if (v == 0) fatal_err("cannot allocate bond vectors", -1);

	blist_temp = (int ***) calloc(num_beads, sizeof(int **));
	if (blist_temp == 0) fatal_err("cannot allocate blist", -1);

	blist_rbc = (int ***) calloc(MRBC, sizeof(int **));
	if (blist_rbc == 0) fatal_err("cannot allocate blist", -1);

	for (n=0 ; n<Nsphere ; n++) 
	{
		centers[n] = (double *) calloc(DIMS, sizeof(double));
		if (centers[n] == 0)  fatal_err ("cannot allocate centers", -1);
	}

	for (i=0 ; i<num_beads ; i++) 
	{
		midpoints[i] = (int *) calloc(MAX_BOND+1, sizeof(int));
		if (midpoints[i] == 0) fatal_err("cannot allocate midpoints", -1);

		blist_temp[i] = (int **) calloc(MAX_BOND+1, sizeof(int *));
		if (blist_temp[i] == 0) fatal_err("cannot allocate blist", -1);

		blist_rbc[i] = (int **) calloc(MAX_BOND+1, sizeof(int *));
		if (blist_rbc[i] == 0) fatal_err("cannot allocate blist", -1);

		for (j=0 ; j<=MAX_BOND ; j++) 
		{
			blist_temp[i][j] = (int *) calloc(DIMS, sizeof(int));
			if (blist_temp[i][j] == 0) fatal_err("cannot allocate blist", -1);

			blist_rbc[i][j] = (int *) calloc(DIMS, sizeof(int));
			if (blist_rbc[i][j] == 0) fatal_err("cannot allocate blist", -1);
		}
	}

	maxNbeads = max(Nbeads[0], Nbeads[1]);

	for(type = 0; type < NTYPES; type++) 
	{
		v[type]= (double **) calloc(maxNbeads, sizeof(double *));
		if (v[type] == 0) fatal_err("cannot allocate bond vectors", -1);

		for(i=0; i<maxNbeads; i++) 
		{
			v[type][i] = (double *) calloc(DIMS, sizeof(double));
			if (v[type][i] == 0) fatal_err("cannot allocate bond vectors", -1);
		}
	}

	if (test)  
		printf("Adding %d spheres, %d particles in total\n", Nsphere, num_beads);

	// Initialize monomer properties
	for(i=0; i < num_beads; i++) 
	{
		if(sphere_pm->ev_type == 0)        /* HS */
			//monomers[i].radius = 0.;
monomers[i].radius=0.5;  // modification 20170313
		else if(sphere_pm->ev_type == 1)   /* WCA */
			monomers[i].radius = 0.5;
		else if(sphere_pm->ev_type == 2)   /* gaussian */
			monomers[i].radius = sphere_pm->Ss;
    else if(sphere_pm->ev_type == 3)
      monomers[i].radius = 0.5;
		else 
		{
			printf("wrong EV type!!\n");
			exit(1);
		}

		monomers[i].rho = 1.0;
		monomers[i].dr2 = 0.0;
		monomers[i].blist[0][0] = 0;
		monomers[i].blist[0][1] = 0;
		monomers[i].blist[0][2] = 0;

		for(j=0; j<DIMS; j++) 
		{
			monomers[i].vel[j]=0.0;
			monomers[i].force[j]=0.0;
			monomers[i].force0[j]=0.0;
			monomers[i].fluid_vel[j]=0.0;
		}
	}

	// Define rest volume and surface area
	for(i=0; i < NTYPES; i++) 
	{
		if(nlevel[i] == -1 || nlevel[i] >=2) 
		{
			if(nlevel[i] == -1)  Vfactor[i] = 0.57;
			else                 Vfactor[i] = 1.0;

			if(sphere_pm->spring == 0) 
			{
				if(sphere_pm->N_per_sphere[i] == 162) 
				{
					var = log(sphere_pm->H_fene[i])/3.0;
					sphere_pm->V0[i] = (59.6+231.0*exp(-0.29*exp(var))) * Vfactor[i];
					sphere_pm->A0[i] = (75.4+138.328*exp(-0.25*exp(var)));
				}
				else if(sphere_pm->N_per_sphere[i] == 642) 
				{
					sphere_pm->V0[i]=(1811.5-187.77*log(sphere_pm->H_fene[i])) * Vfactor[i];
					sphere_pm->A0[i]=(724.14-56.592*log(sphere_pm->H_fene[i]));
				}
			}
			else if(sphere_pm->spring == 2) 
			{
				if(sphere_pm->N_per_sphere[i] == 162) 
				{
					var = log(sphere_pm->H_fene[i])/3.0;
					//sphere_pm->V0[i] = (151.0+155.55*exp(-0.384*exp(var))) * Vfactor[i];
					//sphere_pm->V0[i] = 156.7182;
					//sphere_pm->A0[i] = 141.1168;
          sphere_pm->V0[i] = 153.3733;
          sphere_pm->A0[i] = 138.5641;
          sphere_pm->V0_init[i] = 153.3733;  
					//sphere_pm->A0[i] = (139.6+83.85*exp(-0.388*exp(var)));
          //sphere_pm->V0[i] = 90.51;
					//sphere_pm->A0[i] = 137.57; //140   
				}
				else if(sphere_pm->N_per_sphere[i] == 642) 
				{
					//sphere_pm->V0[i]= (1885.3-102.73*log(sphere_pm->H_fene[i]))*Vfactor[i];
					//sphere_pm->A0[i] = (738.43 - 27.447*log(sphere_pm->H_fene[i]));  
          sphere_pm->V0[i] = 1.2270e+03; // 20170131
					sphere_pm->A0[i] = 5.542563e+02;
					sphere_pm->V0_init[i] = 1.2270e+03;
				}
			}
			else if(sphere_pm->spring == 3) 
			{
				if(sphere_pm->N_per_sphere[i] == 162) 
				{
					if(sphere_pm->H_fene[i] < 200) 
					{
						sphere_pm->V0[i]=(266.45-16.611*log(sphere_pm->H_fene[i]))*Vfactor[i];
						sphere_pm->A0[i] = (201.31-8.8242*log(sphere_pm->H_fene[i]));
					}
					else 
					{
						sphere_pm->V0[i]=(229.84-10.078*log(sphere_pm->H_fene[i]))*Vfactor[i];
						sphere_pm->A0[i] = (179.89 - 5.2169*log(sphere_pm->H_fene[i]));
					}
				}
				else if(sphere_pm->N_per_sphere[i] == 642) 
				{
					sphere_pm->V0[i] = 2100.0 *Vfactor[i];
					sphere_pm->A0[i] = 793.862;
				}
			}
		}
		else if (nlevel[i] == 0)
			sphere_pm->V0[i] = 2.07121;
		else if (nlevel[i] == 1) 
			sphere_pm->V0[i] = 17.1369;
		else 
		{
			fprintf(stderr, "rest volume for nlevel=%d is unknown!", nlevel);
			exit(1);
		}
	}

	// Generate a new config
	/* 1) determine bond length (stretched length) */
	/* 2) start to grow spheres */
	/* 3) write config to a file */
	if(sphere_pm->initconfig == 1 || sphere_pm->initconfig == 2) 
	{   
		// Determine bond length (stretched length) what's the bondlength??
		if(sphere_pm->spring == 0)  /* FENE */
		{
			// bondlength = sphere_pm->Q_fene*(monomers[0].radius*2.0)*MAX_EXT;
			bondlength = 0.98; /* no stretch */
		}
		else if(sphere_pm->spring == 1)  /* WLC */
			bondlength = sphere_pm->nks * sphere_pm->sigma_k * MAX_EXT;
		else if(sphere_pm->spring == 2)  /* harmonic */
//			bondlength = 1.0;  //20170220
      bondlength = 1.6549;
		else {
			printf("wrong spring type!!\n");
			exit(1);
		}
    // Set up the radius of each type of particles
    // what's the relation between radius and bondlength??
		for(i=0; i < NTYPES ; i++) 
		{  
			if(nlevel[i] == -1)
				//radius[i] = 3.30*bondlength;  /* measured relaxation radius */
        radius[i] = 3.91; // 20170307 
			else if (nlevel[i] == 0)
				/* radius of icosahedron (around 0.95*bondlength) */
				radius[i] = bondlength / 4. * sqrt(10.+2.*sqrt(5.)); 
			else if (nlevel[i] == 1) 
				radius[i] = 1.71*bondlength;  /* measured relaxation radius */
			else if (nlevel[i] == 2)
//				radius[i] = 3.30*bondlength;  /* measured relaxation radius */
//        radius[i] = 3.3206*bondlength; // 20170201
        radius[i] = 3.91;
			else if (nlevel[i] == 3)
//				radius[i] = 5.70*bondlength;  /* speculated radius */
        radius[i] = 6.6413*bondlength; // 20170131
			else
				radius[i] = 5.70 * pow(1.98,nlevel[i]-3) * bondlength;
		}

		if(test) 
		  printf("bondlength=%lf,radius[0]=%lf,radius[1]=%lf\n",bondlength,radius[0],radius[1]);
		
		// 2) start to grow spheres                                               
		//    a) set up vertices positions for RBC and sphere                     
		//    b) find a center with enough space around                           
		//    c) randomly rotate the particle, make sure it stays inside the box  

    m = 0; // monomer index
		for(n=0; n < Nsphere; n++)  // Nsphere: the # of spheres
		{     
			// Set up vertices positions for RBC and sphere 
   		if(trial_count > trial_count_max)
			{
				n     = 0;      
				m     = 0;      
				start = 0;      
				type  = 0;      
        radius[0] = radius[0]*0.98;
        radius[1] = radius[1]*0.98;
				if(nlevel[type] == -1)  //i -> type  
					setup_template(v[type], blist_rbc, MRBC, radius[type], work_dir);
				else 
					setup(v, blist_ico, radius);
        /* modification 20170314: 
           'setup_template' is modified so that the radius doesn't affect the particle
           size.
        */
			}
			else
			{
				type = (n < Ntype[0] ? 0 : 1);  // get information about the particle type
				                                // Ntype[]: particle number of each type 
				start = (type==0 ? n * Nbeads[0] : Ntype[0]*Nbeads[0] + (n-Ntype[0])*Nbeads[1]);

				if(nlevel[type] == -1) // Determine vertices position from template file for RBC 
				  setup_template(v[type], blist_rbc, Nbeads[type], radius[type], work_dir);

				else // Determine vertices position and bonding list for the icosahedron
				  setup(v, blist_ico, radius);				
			}

			if(test) { printf("Adding sphere %d\n", n); }

			// Find a center with enough space around  
			trial_count = 0;
			do 
			{
				if(trial_count > trial_count_max) 
				{
					if(n_proc == 0) 
					{
						printf("Tried %d times but cannnot insert sphere %d!\n", trial_count_max, n);
						break;  // Leave the do-while loop and go to the if-statement
					}         // below, finally come back to the for-loop above.  
				}

				// Choose center of the sphere randomly between 0 and maxsize[j]-1.
        /* modification 20170314:
           In the particle coordinate system, the ends of the box is 0 and max-1. 
           This statement mightbe WRONG!!!

           modification 20170316: 
        */
				for(j=0; j < DIMS; j++)
				{
// modified 20160915
  			  // centers[n][j] = ( (double)rand()/(double)RAND_MAX ) * (maxsize[j]-2.0) + 1.;
					//(double)rand() / ( (double)(RAND_MAX) + 1.0 ) * (maxsize[j]-2.0) + 1.0;

					centers[n][j] = ((double)rand()/(double)RAND_MAX) * (maxsize[j]); //20170316
          //centers[n][j] = ((double)rand()/(double)RAND_MAX) * (maxsize[j]-2*(1.95+0.5)) + (1.95 + 0.5);
          //centers[n][j] = ((double)rand()/(double)RAND_MAX) * (maxsize[j]-2*(1.95+0.25)) + (1.95 + 0.25);
          //double position = 22 + (rand() % ( (int)(10*(maxsize[j]-2.2)) - 22 + 1));
          //centers[n][j] = position/10.;
				}
        /*modification 20170314:
          This setup should match the wall force function. The particle radius is 
          reduced by half, so the factor 3.15 = 3.91*0.5 + factor, where factor = 1.2 
          taken from the wall force function. The upper limit is max-1-1.2-3.91/2 = 
          max - 4.15 
        */
        if(wall_flag==1)
 				  centers[n][1] = ((double)rand()/(double)RAND_MAX) * ((maxsize[1]-1)-2*(1.2+1.95+0.)) + 
                          1.2+1.95+0.;

        /* modification 20170314: 
           the fluid domain goes from 1 to max. The corresponding particle coordinate 
           is from 0 to max-1. To locate the particle in the middle of the fluid 
           domain, the coordinate should be (max-1)/2 in the particle coordinate system
        */
				if(sphere_pm->initconfig == 2 && Nsphere == 1) 
					for(j=0; j<DIMS; j++)
						//centers[n][j]=(maxsize[j] - 1.0 - 2.0) / 2.0 + 1.0;
            centers[n][j] = (maxsize[j] - 1.0) / 2.0;

        if(sphere_pm->initconfig == 2 && Nsphere == 2) 
        {
          printf("Temporary design! line 381 in init_sphere.c should be modified!\n");
          for(j=0; j<DIMS; j++)
					  centers[n][j]=(maxsize[j] - 1.0 - 2.0) / 2.0 + 1.0;
          centers[0][1]+=3;
          centers[1][1]-=3;
          centers[0][0]=5;
          centers[1][0]+=15;
        }
        // 20161220

				// Check sphere-sphere overlap 
				overlap = FALSE;
				trial_count++;
				for(n1=0 ; n1 < n ; n1++) 
				{
					type1 = (n1 < Ntype[0] ? 0 : 1);  
					dx = centers[n][0] - centers[n1][0];
					dy = centers[n][1] - centers[n1][1];
					dz = centers[n][2] - centers[n1][2];
					if (wall_flag < 3)        /* no x-wall */
						dx = n_image(dx,max_x);
					if (wall_flag < 2)        /* no z-wall */
						dz = n_image(dz,max_z);
					if (wall_flag < 1)        /* no y-wall */
						dy = n_image(dy,max_y);
					r2 = dx*dx + dy*dy + dz*dz;
          /* modification 20170314: 
             Now the initilization process is only suitable for particles with a radius 
             of 3.91. At the beginning of the initialization process, the radius is 
             reduced by half. 4.96 = (3.91/2) + (3.91/2) + factor, where factor=1.05;
          */
					if(sqrt(r2) < 4.96) { 
					//if(sqrt(r2) <= (radius[type] + radius[type1]) * 1.0 + 1.06) {      // Why 1.06?
				    overlap = TRUE;
            overlap_pp++;
					}
				}

				// Check wall-overlap  // Why 1.05???
				if(wall_flag >= 3)  /* x-wall */
				{  
					if(centers[n][0] + radius[type]*1.0 + 1.05 > max_x - 1 || 
						 centers[n][0] - radius[type]*1.0 - 1.05 < 0) {  
						overlap = TRUE;
					}	  
				}
				if(wall_flag >= 2)  /* z-wall */
				{  
					if(centers[n][2] + radius[type]*1.0 + 1.05 > max_z - 1 || 
						 centers[n][2] - radius[type]*1.0 - 1.05 < 0) {  
						overlap = TRUE;
					}
				}
				if(wall_flag >= 1)  /* y-wall */
				{ 
          /* modification 20170314:
             This criterion should match the setup of the wall force. The particle 
             radius is reduced, and I use the reduced radius here rather than the 
             expected one because I want to save space. I am not sure if this setup 
             leads to any instability problem due to the wall force acting on 
             particles.The factor 1.2 is taken from the wall force function.
          */ 
				  //if(centers[n][1] + radius[type]*1.0 + 1.05 > max_y - 1 || centers[n][1] - radius[type]*1.0 - 1.05 < 0) {
				  if((centers[n][1] + 1.95) > ((max_y-1)-1.2-0.) || (centers[n][1]-1.95) < (1.2+0.)) {
						overlap = TRUE;
            overlap_pw++;
					}
				}
			}while(overlap == TRUE);

			if(trial_count > trial_count_max) {
        n=0;                           // Modification 20170316!!!
        continue;  // Come back to the for-loop above.
      }

			if(test)
				printf("sphere %d centers at (%le,%le,%le)\n", n, centers[n][0], centers[n][1], 
        centers[n][2]);

			// Randomly rotate the particle, make sure it stays inside the box

			//rot_flag = 1;        Modification 20170320
			//while(rot_flag == 1)
      do 
			{
        rot_flag = 0; // Modification 20170320
				/* choose a direction for the z-axis of the sphere */
				phi = (double)rand() / (double)(RAND_MAX) * 2. * M_PI;
				theta = (double)rand() / (double)(RAND_MAX) * M_PI;
//       theta= M_PI * 90. / 180;
//       phi = 0.5*M_PI;
				/* z axis */
				axis[2][0] = cos(phi)*sin(theta);
				axis[2][1] = sin(phi)*sin(theta);
				axis[2][2] = cos(theta);
				/* y axis is chosen arbitrarily by rotating z-axis down with pi/2 */
				axis[1][0] = cos(phi)*sin(theta+M_PI/2.);
				axis[1][1] = sin(phi)*sin(theta+M_PI/2.);
				axis[1][2] = cos(theta+M_PI/2.);
				/* x axis is dicided by a cross product */
				product(axis[1], axis[2], axis[0]);
//       axis[0][0]=cos(theta);  axis[0][1]=-sin(theta); axis[0][2]=0.;
//       axis[1][0]=sin(theta);  axis[1][1]=cos(theta);  axis[1][2]=0.;
//       axis[2][0]=0.;          axis[2][1]=0.;          axis[2][2]=1.;

				/* normalize again in case there is numerical error */
				for(j=0; j < DIMS; j++) 
        {
					r2 = 0.;
					for(k=0; k < DIMS; k++)
						r2 += axis[j][k]*axis[j][k];
					r2 = sqrt(r2);
					for(k=0; k < DIMS; k++)
						axis[j][k] /= r2;
				}
				/* insert a template */
				if(nlevel[type] == -1) 
				{
					for(i=0; i < Nbeads[type]; i++) 
					{
            // Modification 20170320
            if(rot_flag==1)
              continue;
						m = start + i;
						for(j=0; j < DIMS; j++) 
						{
							monomers[m].pos_pbc[j] = centers[n][j];
							for(k=0; k < DIMS; k++) 
						    monomers[m].pos_pbc[j] += v[type][i][k]*axis[k][j];               
							monomers[m].pos[j] = monomers[m].pos_pbc[j];
						}
						monomers[m].pos[0] = box(monomers[m].pos_pbc[0],maxsize[0]);
						if(wall_flag < 2) 
							monomers[m].pos[2] = box(monomers[m].pos_pbc[2],maxsize[2]);
						if(wall_flag < 1)
							monomers[m].pos[1] = box(monomers[m].pos_pbc[1],maxsize[1]);

						/* check if the rotated vertex position go out of walls */
						if(wall_flag == 0)
						  rot_flag = 0;

						for(j=1; j<=2; j++) 
						{
							if(wall_flag == j) 
							{
								//for(i=1; i<=j; i++) // it will overwrite i 
								for(int dim=1; dim <= j; dim++) 
								{
                  /* modification 20170316:

                  */
									if(monomers[m].pos_pbc[dim] > (maxsize[dim]-1.0-1.2-1.) || 
										 monomers[m].pos_pbc[dim] < (0+1.2+1.)) 
									{
										rot_flag = 1;
                    rotation++;
										break;
									}
									else 
										rot_flag = 0;
								}
							}
						}
            if(rot_flag==1)  // Modification 20170320
              continue;

						monomers[m].blist[0][0] = blist_rbc[i][0][0];

						for (j=1 ; j<=blist_rbc[i][0][0] ; j++) 
						{
							monomers[m].blist[j][0] = start + blist_rbc[i][j][0];
							monomers[m].blist[j][1] = start + blist_rbc[i][j][1];
							monomers[m].blist[j][2] = start + blist_rbc[i][j][2];
						}

						monomers[m].sphere_id = n;
					}
				}
				// insert an icosahedron 
        // Modification 20170320: I modified the rotation processs for RBC. I am not sure if it also
        // works for spheres.
				else 
				{
					m = start;   
					rot_flag = 0;  // for a sphere, we only need to rotate monomer positions once,
          // so rot_flag is set to be zero.
					for(i=0; i < NVERTICES; i++) 
					{
            // rotate monomer positions
						for(j=0; j < DIMS; j++) { 
							monomers[m].pos_pbc[j] = centers[n][j];
							for(k=0; k < DIMS; k++)
								monomers[m].pos_pbc[j] += v[type][i][k] * axis[k][j];
							monomers[m].pos[j] = box(monomers[m].pos_pbc[j],maxsize[j]);
						}

						bcount = 0;
						for(j=1; j <= NBONDS; j++) 
						{
							k = start + blist_ico[i][j][0];  // m is bonded to k
							// Add k to the blist of m only if m is not already in the blist of k
              // It's clever design for checking the bonded list. If I exchange the order 
              // of k and m, i.e. check(m,k,monomers), I can't correctly check the list.
							if(!check_bond(k,m,monomers)) // If m is NOT already in the blist of k
							{
								bcount++;
								monomers[m].blist[bcount][0] = k;
								monomers[m].blist[bcount][1] = start + blist_ico[i][j][1];
								monomers[m].blist[bcount][2] = start + blist_ico[i][j][2];
							}
						}
						monomers[m].blist[0][0] = bcount;
						monomers[m].sphere_id = n;
						monomers[m].type = type;
						m++;
					}
					// sphere triangulation  
					for(level=0; level < nlevel[type]; level++) 
					{
						m_temp = m;

						// refine the mesh by creating additional points 
						for(n1=start; n1 < m_temp; n1++)
						{
							for(j=1; j <= monomers[n1].blist[0][0]; j++) // loop over all bonds
							{  
								n2 = monomers[n1].blist[j][0];

								// locate the midpoint 
								r2 = 0.;
								for(k=0; k<DIMS; k++) 
								{
									p_mid[k] = (monomers[n1].pos_pbc[k] + monomers[n2].pos_pbc[k])/2.;
									dr[k] = p_mid[k] - centers[n][k];
									r2 += dr[k]*dr[k];
								}
								r2 = sqrt(r2);

								// push out the midpoint to the sphere 
								for (k=0 ; k < DIMS ; k++) 
								{
									monomers[m].pos_pbc[k] = centers[n][k] + dr[k]*radius[type]/r2;
									monomers[m].pos[k] = box(monomers[m].pos_pbc[k], maxsize[k]);
								}
								monomers[m].sphere_id = n;
								monomers[m].type = type;
								midpoints[n1][j] = m;
								//printf("n=%d start=%d m=%d monpos=%le %le %le\n", 
								//n, start, m, monomers[m].pos_pbc[0],
								//monomers[m].pos_pbc[1],monomers[m].pos_pbc[2]);
								m++;
							}
						}
            
            // set up blist for the new points and update blist for the existing points
            // This block is very sophisticated. I don't think it's easy to write a 
            // better version by myself. // 20170218
						for(n1 = start; n1 < m_temp; n1++) {
							for(j=1; j <= monomers[n1].blist[0][0]; j++) 
							{ 
								// This is the relative position of the points (facing outside) 

								n2 = monomers[n1].blist[j][0];                    /*              4               */
								n3 = monomers[n1].blist[j][1];                    /*              /\              */
								n4 = monomers[n1].blist[j][2];                    /*          14 /__\ 24          */
                                                  								/*            /\  /\            */
								n12 = midpoints[n1][j];                           /*        1  /__12__\  2        */
								n13 = get_midpoint(n1, n3, monomers, midpoints);  /*           \  /\  /           */
								n14 = get_midpoint(n1, n4, monomers, midpoints);  /*            \/__\/ 23         */
								n23 = get_midpoint(n2, n3, monomers, midpoints);  /*          13 \  /             */
								n24 = get_midpoint(n2, n4, monomers, midpoints);  /*              \/              */
								                                                  /*              3               */

								/* update blist: the n1-n2 bond now becomes n1-n12 bond */
								/* we still need the old blist for getting midpoints. 
									 So we store the update in blist_temp and recover it later */
								blist_temp[n1][j][0] = n12;
								blist_temp[n1][j][1] = n13;
								blist_temp[n1][j][2] = n14;

								/* set up blist for n12 */
								bcount = 1;
								monomers[n12].blist[bcount][0] = n2;
								monomers[n12].blist[bcount][1] = n23;
								monomers[n12].blist[bcount][2] = n24;
								if (!check_bond(n13,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n13;
									monomers[n12].blist[bcount][1] = n1;
									monomers[n12].blist[bcount][2] = n23;
								}
								if (!check_bond(n23,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n23;
									monomers[n12].blist[bcount][1] = n13;
									monomers[n12].blist[bcount][2] = n2;
								}
								if (!check_bond(n14,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n14;
									monomers[n12].blist[bcount][1] = n24;
									monomers[n12].blist[bcount][2] = n1;
								}
								if (!check_bond(n24,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n24;
									monomers[n12].blist[bcount][1] = n2;
									monomers[n12].blist[bcount][2] = n14;
								}
								monomers[n12].blist[0][0] = bcount;
							}
            }
						/* recover blist for the old points */
						for (n1=start ; n1<m_temp ; n1++)
							for (j=1 ; j<=monomers[n1].blist[0][0] ; j++)
								for (k=0 ; k<3 ; k++)
									monomers[n1].blist[j][k] = blist_temp[n1][j][k];
					}
				}
			}while(rot_flag==1); // Modification 20170320
		}
		// Write config to a file
		if(test) 
    {
		  file_name(test_file, work_dir, -1);
			file_ptr = fopen(test_file, "w");
			if(file_ptr == 0)  
			  fatal_err("Could not open test.dat", -1);
			for(m=0; m < num_beads; m++) {
        fprintf(file_ptr,"%le %le %le\n",monomers[m].pos[0],monomers[m].pos[1], 
                monomers[m].pos[2]);	
      }
			fclose(file_ptr);
		}
    // Modification 20170320: Record the initialization process.
    sprintf(filename,"%s/data/initializationRecord.dat",work_dir); 
    stream = fopen(filename,"w");  
    fprintf(stream,"particle-particle overlap: %d\n", overlap_pp);
    fprintf(stream,"particle-wall overlap:     %d\n", overlap_pw);
    fprintf(stream,"rotation:                  %d\n", rotation);
    fclose(stream);
	}
	else if (sphere_pm->initconfig == 3)
	{
		/////////////////////////////////////////////////
		//                                             //
		// 1) determine bond length (stretched length) //
		//                                             //
		/////////////////////////////////////////////////

		if(sphere_pm->spring == 0)  /* FENE */
		{
			// bondlength = sphere_pm->Q_fene*(monomers[0].radius*2.0)*MAX_EXT;
			bondlength = 0.98; /* no stretch */
		}
		else if(sphere_pm->spring == 1)  /* WLC */
			bondlength = sphere_pm->nks * sphere_pm->sigma_k * MAX_EXT;

		else if(sphere_pm->spring == 2)  /* harmonic */
			bondlength = 1.0;

		else 
		{
			printf("wrong spring type!!\n");
			exit(1);
		}

		for (i=0; i < NTYPES ; i++) 
		{
			if (nlevel[i] == -1)
				radius[i] = 3.30*bondlength;  /* measured relaxation radius */

			else if (nlevel[i] == 0)
				/* radius of icosahedron (around 0.95*bondlength) */
				radius[i] = bondlength / 4. * sqrt(10.+2.*sqrt(5.)); 

			else if (nlevel[i] == 1) 
				radius[i] = 1.71*bondlength;  /* measured relaxation radius */

			else if (nlevel[i] == 2)
				radius[i] = 3.30*bondlength;  /* measured relaxation radius */

			else if (nlevel[i] == 3)
				radius[i] = 5.70*bondlength;  /* speculated radius */

			else
				radius[i] = 5.70 * pow(1.98,nlevel[i]-3) * bondlength;
		}

		if (test)
		{
			printf("bondlength=%lf, radius[0]=%lf, radius[1]=%lf\n", 
					bondlength, radius[0], radius[1]);
		}

		/* monomer index */
		m = 0;
		char filedir[500];
		sprintf(filedir, "%s/init/HCP_config.dat", work_dir);
		FILE *stream_HCP;
		stream_HCP = fopen(filedir,"r");   

		////////////////////////////////////////////////////////////////////////////
		//                                                                        //
		// 2) start to grow spheres                                               //
		//    a) set up vertices positions for RBC and sphere                     //
		//    b) find a center with enough space around                           //
		//    c) randomly rotate the particle, make sure it stays inside the box  //
		//                                                                        //
		////////////////////////////////////////////////////////////////////////////

		for(n=0; n < Nsphere; n++)    // Nsphere: total particle number 
		{     
			if (test)    {printf("Adding sphere %d\n", n);}

			/* a. set up vertices positions for RBC and sphere ----------------------------*/

			type = (n < Ntype[0] ? 0 : 1);  // get information about the particle type
			// Ntype[]: particle number of each type 
			start = (type==0 ? 
					n * Nbeads[0] : Ntype[0]*Nbeads[0] + (n-Ntype[0])*Nbeads[1]);

			if (overlap == TRUE)
			{
				radius[0] = 0.98*radius[0];
				radius[1] = 0.98*radius[1];

				if(nlevel[type] == -1){
					/* determine vertices position from template file for RBC */
					setup_template(v[type], blist_rbc, Nbeads[type], radius[type], work_dir);
				}
				else{
					/* determine vertices position and bonding list for the icosahedron */
					setup(v, blist_ico, radius);
				}
			}
			else
			{
				if(nlevel[type] == -1){
					/* determine vertices position from template file for RBC */
					setup_template(v[type], blist_rbc, Nbeads[type], radius[type], work_dir);
				}
				else{
					/* determine vertices position and bonding list for the icosahedron */
					setup(v, blist_ico, radius);
				}
			}

			/*  b. find a center with enough space around --------------------------------*/ 

			fscanf(stream_HCP, "%le %le %le", 
					&centers[n][0], &centers[n][1], &centers[n][2]);

			//-----------------------------
			// check sphere-sphere overlap 
			//-----------------------------

			overlap = FALSE;

			for (n1=0 ; n1 < n ; n1++) 
			{
				type1 = (n1 < Ntype[0] ? 0 : 1);  

				dx = centers[n][0] - centers[n1][0];
				dy = centers[n][1] - centers[n1][1];
				dz = centers[n][2] - centers[n1][2];

				if (wall_flag < 3)        /* no x-wall */
					dx = n_image(dx,max_x);
				if (wall_flag < 2)        /* no z-wall */
					dz = n_image(dz,max_z);
				if (wall_flag < 1)        /* no y-wall */
					dy = n_image(dy,max_y);

				r2 = dx*dx + dy*dy + dz*dz;

				if (sqrt(r2) <= (radius[type] + radius[type1]) * 1.0 + 1.06)
				{  
					overlap = TRUE;
					printf("sphere-sphere overlap\n");
					//PAUSE;
				}
			}

			//-----------------------------------
			// check wall-overlap  // Why 1.05???
			//-----------------------------------

			if (wall_flag >= 3)  /* x-wall */
				if (centers[n][0] + radius[type]*1.0 + 1.05 > max_x - 1 || 
						centers[n][0] - radius[type]*1.0 - 1.05 < 0)  

					overlap = TRUE;

			if (wall_flag >= 2)  /* z-wall */
				if (centers[n][2] + radius[type]*1.0 + 1.05 > max_z - 1 || 
						centers[n][2] - radius[type]*1.0 - 1.05 < 0) 

					overlap = TRUE;

			if (wall_flag >= 1)  /* y-wall */
				if (centers[n][1] + radius[type]*1.0 + 1.05 > max_y - 1 || 
						centers[n][1] - radius[type]*1.0 - 1.05 < 0)
				{  
					overlap = TRUE;
					printf("ywall-sphere overlap\n");
					//printf("upper: center=%lf radius=%lf result=%lf\n", centers[n][1], radius[type], centers[n][1] + radius[type]*1.0 + 1.05);
					//printf("lower: center=%lf radius=%lf result=%lf\n", centers[n][1], radius[type], centers[n][1] - radius[type]*1.0 - 1.05);

					//PAUSE;
				}

			if(overlap == TRUE)  
			{
				n = -1;
				m = 0;
				continue;  // Come back to the for-loop above.
			}

			if (test)
			{
				printf("sphere %d centers at (%le,%le,%le)\n", 
						n, centers[n][0], centers[n][1], centers[n][2]);
			}

			/* c. randomly rotate the particle, make sure it stays inside the box --------*/ 

			rot_flag = 1;
			while(rot_flag == 1) 
			{
				/* choose a direction for the z-axis of the sphere */
				phi = (double)rand() / (double)(RAND_MAX) * 2. * M_PI;
				theta = (double)rand() / (double)(RAND_MAX) * M_PI;

//       theta= M_PI * 90. / 180;
//       phi = 0.5*M_PI;

				/* z axis */
				axis[2][0] = cos(phi)*sin(theta);
				axis[2][1] = sin(phi)*sin(theta);
				axis[2][2] = cos(theta);

				/* y axis is chosen arbitrarily by rotating z-axis down with pi/2 */
				axis[1][0] = cos(phi)*sin(theta+M_PI/2.);
				axis[1][1] = sin(phi)*sin(theta+M_PI/2.);
				axis[1][2] = cos(theta+M_PI/2.);

				/* x axis is dicided by a cross product */
				product(axis[1], axis[2], axis[0]);

//       axis[0][0]=cos(theta);  axis[0][1]=-sin(theta); axis[0][2]=0.;
//       axis[1][0]=sin(theta);  axis[1][1]=cos(theta);  axis[1][2]=0.;
//       axis[2][0]=0.;          axis[2][1]=0.;          axis[2][2]=1.;

				/* normalize again in case there is numerical error */
				for (j=0 ; j<DIMS ; j++) 
				{
					r2 = 0.;

					for (k=0 ; k<DIMS ; k++)
						r2 += axis[j][k]*axis[j][k];

					r2 = sqrt(r2);

					for (k=0 ; k<DIMS ; k++)
						axis[j][k] /= r2;
				}

				/* insert a template */

				if(nlevel[type] == -1) 
				{
					for (i=0 ; i < Nbeads[type] ; i++) 
					{
						m = start + i;
						for (j=0 ; j < DIMS ; j++) 
						{
							monomers[m].pos_pbc[j] = centers[n][j];

							for (k=0 ; k<DIMS ; k++) 
								monomers[m].pos_pbc[j] += v[type][i][k]*axis[k][j];

							monomers[m].pos[j] = monomers[m].pos_pbc[j];
						}

						monomers[m].pos[0] = box(monomers[m].pos_pbc[0],maxsize[0]);

						if(wall_flag < 2) 
							monomers[m].pos[2] = box(monomers[m].pos_pbc[2],maxsize[2]);

						if(wall_flag < 1)
							monomers[m].pos[1] = box(monomers[m].pos_pbc[1],maxsize[1]);

						/* check if the rotated vertex position go out of walls */
						if(wall_flag ==0)
							rot_flag = 0;

						for(j=1; j<=2; j++) 
						{
							if(wall_flag == j) 
							{
								for(i=1; i<=j; i++) 
								{
									if (monomers[m].pos_pbc[i] > maxsize[i] - 1.0 || monomers[m].pos_pbc[i] < 0) 
									{
										rot_flag = 1;
										break;
									}
									else 
										rot_flag = 0;
								}
							}
						}

						monomers[m].blist[0][0] = blist_rbc[i][0][0];

						for (j=1 ; j<=blist_rbc[i][0][0] ; j++) 
						{
							monomers[m].blist[j][0] = start + blist_rbc[i][j][0];
							monomers[m].blist[j][1] = start + blist_rbc[i][j][1];
							monomers[m].blist[j][2] = start + blist_rbc[i][j][2];
						}

						monomers[m].sphere_id = n;
					}
				}

				/* insert an icosahedron */
				else 
				{
					m = start;    // ?
					rot_flag = 0;

					for (i=0; i < NVERTICES; i++) 
					{
						for (j=0; j < DIMS; j++) 
						{
							monomers[m].pos_pbc[j] = centers[n][j];

							for (k=0; k < DIMS; k++)
								monomers[m].pos_pbc[j] += v[type][i][k]*axis[k][j];

							monomers[m].pos[j] = box(monomers[m].pos_pbc[j],maxsize[j]);
						}

						bcount = 0;

						for (j=1; j <= NBONDS; j++) 
						{
							k = start + blist_ico[i][j][0];  /* m is bonded to k */
							/* add k to the blist of m only if m is not already in the blist of k */
							if (!check_bond(k,m,monomers)) 
							{
								bcount++;
								monomers[m].blist[bcount][0] = k;
								monomers[m].blist[bcount][1] = start + blist_ico[i][j][1];
								monomers[m].blist[bcount][2] = start + blist_ico[i][j][2];
							}
						}

						monomers[m].blist[0][0] = bcount;
						monomers[m].sphere_id = n;
						monomers[m].type = type;
						// printf("n=%d start=%d m=%d monpos=%le %le %le\n", 
						//        n, start, m, monomers[m].pos_pbc[0], monomers[m].pos_pbc[1], 
						//        monomers[m].pos_pbc[2]);
						m++;
					}

					//--------------------------------------------------  
					// sphere triangulation  
					//--------------------------------------------------

					for (level=0; level < nlevel[type]; level++) 
					{
						m_temp = m;

						//----------------------------------------------
						// refine the mesh by creating additional points 
						//----------------------------------------------

						for (n1=start; n1 < m_temp; n1++)
						{
							for (j=1; j <= monomers[n1].blist[0][0]; j++) /* loop over all bonds */
							{  
								n2 = monomers[n1].blist[j][0];

								//--------------------
								// locate the midpoint 
								//--------------------		

								r2 = 0.;
								for (k=0; k<DIMS; k++) 
								{
									p_mid[k] = (monomers[n1].pos_pbc[k]+monomers[n2].pos_pbc[k])/2.;
									dr[k] = p_mid[k] - centers[n][k];
									r2 += dr[k]*dr[k];
								}
								r2 = sqrt(r2);

								//------------------------------------
								// push out the midpoint to the sphere 
								//------------------------------------

								for (k=0 ; k < DIMS ; k++) 
								{
									monomers[m].pos_pbc[k] = centers[n][k] + dr[k]*radius[type]/r2;
									monomers[m].pos[k] = box(monomers[m].pos_pbc[k], maxsize[k]);
								}
								monomers[m].sphere_id = n;
								monomers[m].type = type;
								midpoints[n1][j] = m;
								//printf("n=%d start=%d m=%d monpos=%le %le %le\n", 
								//n, start, m, monomers[m].pos_pbc[0],
								//monomers[m].pos_pbc[1],monomers[m].pos_pbc[2]);
								m++;
							}
						}

						//-------------------------------------------------------------------------
						// set up blist for the new points and update blist for the existing points 
						//-------------------------------------------------------------------------

						for (n1=start ; n1<m_temp ; n1++)
							for (j=1 ; j<=monomers[n1].blist[0][0] ; j++) 
							{  /* loop over all bonds */

								/* This is the relative position of the points (facing outside) */

								n2 = monomers[n1].blist[j][0];                    /*              4               */
								n3 = monomers[n1].blist[j][1];                    /*              /\              */
								n4 = monomers[n1].blist[j][2];                    /*          14 /__\ 24          */
								/*            /\  /\            */
								n12 = midpoints[n1][j];                           /*        1  /__12__\  2        */
								n13 = get_midpoint(n1, n3, monomers, midpoints);  /*           \  /\  /           */
								n14 = get_midpoint(n1, n4, monomers, midpoints);  /*            \/__\/ 23         */
								n23 = get_midpoint(n2, n3, monomers, midpoints);  /*          13 \  /             */
								n24 = get_midpoint(n2, n4, monomers, midpoints);  /*              \/              */
								/*              3               */
								/*                              */

								/* update blist: the n1-n2 bond now becomes n1-n12 bond */
								/* we still need the old blist for getting midpoints. 
									 So we store the update in blist_temp and recover it later */
								blist_temp[n1][j][0] = n12;
								blist_temp[n1][j][1] = n13;
								blist_temp[n1][j][2] = n14;

								/* set up blist for n12 */
								bcount = 1;
								monomers[n12].blist[bcount][0] = n2;
								monomers[n12].blist[bcount][1] = n23;
								monomers[n12].blist[bcount][2] = n24;
								if (!check_bond(n13,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n13;
									monomers[n12].blist[bcount][1] = n1;
									monomers[n12].blist[bcount][2] = n23;
								}
								if (!check_bond(n23,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n23;
									monomers[n12].blist[bcount][1] = n13;
									monomers[n12].blist[bcount][2] = n2;
								}
								if (!check_bond(n14,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n14;
									monomers[n12].blist[bcount][1] = n24;
									monomers[n12].blist[bcount][2] = n1;
								}
								if (!check_bond(n24,n12,monomers)) {
									bcount++;
									monomers[n12].blist[bcount][0] = n24;
									monomers[n12].blist[bcount][1] = n2;
									monomers[n12].blist[bcount][2] = n14;
								}
								monomers[n12].blist[0][0] = bcount;
							}

						/* recover blist for the old points */
						for (n1=start ; n1<m_temp ; n1++)
							for (j=1 ; j<=monomers[n1].blist[0][0] ; j++)
								for (k=0 ; k<3 ; k++)
									monomers[n1].blist[j][k] = blist_temp[n1][j][k];
					}
				}
			}
		}

		/* 3) write config to a file */
		if (test) 
		{
			file_name (test_file, work_dir, -1);

			file_ptr = fopen (test_file, "w");

			if (file_ptr == 0)  
				fatal_err("Could not open test.dat", -1);

			for (m=0 ; m<num_beads ; m++)
			{
				fprintf(file_ptr, "%le %le %le\n", 
						monomers[m].pos[0], monomers[m].pos[1], monomers[m].pos[2]);
			}

			fclose(file_ptr);
		}
	}

	/* read in the initial sphere configuration from file -----------------------*/

	else if(sphere_pm->initconfig == 4) 
	{
		sprintf(filename, "%s/init/init.config", work_dir);
		printf("init config file %s\n", filename);
		stream = fopen(filename, "r");      
		fscanf(stream, "%d %d %d %d %d", &sphere_pm->Ntype[0], &sphere_pm->Ntype[1], 
				&sphere_pm->nlevel[0], &sphere_pm->nlevel[1], &sphere_pm->num_beads);

		sphere_pm->Nsphere = sphere_pm->Ntype[0] + sphere_pm->Ntype[1];
		Nsphere = sphere_pm->Nsphere;
		num_beads = sphere_pm->num_beads;

		for (j=0 ; j<NTYPES ; j++) 
		{
			sphere_pm->N_per_sphere[j] = 12;
			sphere_pm->face_per_sphere[j] = 20;
			temp = 30;
			for (i=0 ; i<sphere_pm->nlevel[j] ; i++) 
			{
				sphere_pm->N_per_sphere[j] += temp;
				sphere_pm->face_per_sphere[j] *= 4;
				temp *= 4;
			}
		}

		for (i=0 ; i<NTYPES ; i++) 
		{
			Ntype[i] = sphere_pm->Ntype[i];
			nlevel[i] = sphere_pm->nlevel[i];
			Nfaces[i] = sphere_pm->face_per_sphere[i];
			Nbeads[i] = sphere_pm->N_per_sphere[i];
		}

		for(i=0; i<Nsphere; i++) 
			fscanf(stream, "%le %le", &spheres[i].disp2, &spheres[i].dr2);

		for(i=0; i<Nsphere; i++) 
		{
			type = (i<Ntype[0] ? 0 : 1);
			start = (type == 0 ? i*Nbeads[0] : Ntype[0]*Nbeads[0] + (i-Ntype[0])*Nbeads[1]);
			if(n_proc == 0)
				printf("#sphere %d \n", i);
			for(n=0; n<Nbeads[type]; n++) 
			{
				pnum = start+n;
				fscanf(stream, "%le %le %le %le %le %le %le %le %le %le", 
						&monomers[pnum].pos[0], &monomers[pnum].pos[1], 
						&monomers[pnum].pos[2], &monomers[pnum].pos_pbc[0], 
						&monomers[pnum].pos_pbc[1], &monomers[pnum].pos_pbc[2], 
						&monomers[pnum].vel[0], &monomers[pnum].vel[1], 
						&monomers[pnum].vel[2], &monomers[pnum].dr2);

				monomers[pnum].sphere_id=i;
				monomers[pnum].type = type;

				if(n_proc == 0)
					printf("#monomer %d at (%le %le %le) (%le %le %le)\n", 
							pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], 
							monomers[pnum].pos[2], monomers[pnum].pos_pbc[0], 
							monomers[pnum].pos_pbc[1], monomers[pnum].pos_pbc[2]);
			}
		}

		fclose(stream);

		/*read in bond list*/
		sprintf(filename, "%s/init/bond.dat", work_dir);
		stream = fopen(filename, "r");
		for(i=0; i<Nbeads[0]; i++) 
		{
			fscanf(stream, "%*s %d\n",  &monomers[i].blist[0][0]);
			for(j=1; j<=monomers[i].blist[0][0]; j++)
				fscanf(stream, "%d %d %d\n", &monomers[i].blist[j][0], 
						&monomers[i].blist[j][1], &monomers[i].blist[j][2]);
		}
		for(i=0; i<Nbeads[1]; i++) 
		{
			n=Ntype[0]*Nbeads[0]+i;
			fscanf(stream, "%*s %d\n",  &monomers[n].blist[0][0]);
			for(j=1; j<=monomers[n].blist[0][0]; j++)
				fscanf(stream, "%d %d %d\n", &monomers[n].blist[j][0], 
						&monomers[n].blist[j][1], &monomers[n].blist[j][2]);
		}
		fclose(stream);

		for(i=1; i<Nsphere; i++)  
		{
			type = (i<Ntype[0] ? 0 : 1);
			start0 = (type == 0 ? 0 : Ntype[0]*Nbeads[0]);
			start = (type == 0 ? i*Nbeads[0] : Ntype[0]*Nbeads[0] + (i-Ntype[0])*Nbeads[1]);
			for(j=0; j<Nbeads[type]; j++) 
			{
				monomers[start+j].blist[0][0] = monomers[start0+j].blist[0][0];
				for(n=1; n<=monomers[j].blist[0][0]; n++) 
				{
					monomers[start+j].blist[n][0] = start+(monomers[start0+j].blist[n][0]-start0);
					monomers[start+j].blist[n][1] = start+(monomers[start0+j].blist[n][1]-start0);
					monomers[start+j].blist[n][2] = start+(monomers[start0+j].blist[n][2]-start0);
				}
			}
		}
	}

	/*
		 sprintf(filename, "%s/init/bond.dat", work_dir);
		 stream = fopen(filename, "w");
		 for(i=0; i<Nbeads[0]; i++) {
		 fprintf(stream, "%d %d\n",  i, monomers[i].blist[0][0]);
		 for(j=1; j<=monomers[i].blist[0][0]; j++)
		 fprintf(stream, "%d %d %d\n", monomers[i].blist[j][0], 
		 monomers[i].blist[j][1], monomers[i].blist[j][2]);
		 }

		 for(i=0; i<Nbeads[1]; i++) {
		 n=Ntype[0]*Nbeads[0]+i;
		 fprintf(stream, "%d %d\n",  n, monomers[n].blist[0][0]);
		 for(j=1; j<=monomers[n].blist[0][0]; j++)
		 fprintf(stream, "%d %d %d\n", monomers[n].blist[j][0], 
		 monomers[n].blist[j][1], monomers[n].blist[j][2]);
		 }
		 fclose(stream);
	 */

	// add faces 
/*
	f=0;
	for(i=0; i < num_beads; i++) 
	{
		for(j=1; j <= monomers[i].blist[0][0]; j++) 
		{
			if(!check_face(i, monomers[i].blist[j][0], monomers[i].blist[j][1], faces, f)) 
			{
				if(i < Ntype[0]*Nbeads[0])
					faces[f].sphere_id = i/Nbeads[0];
				else
					faces[f].sphere_id = Ntype[0]+(i-Ntype[0]*Nbeads[0])/Nbeads[1];

				faces[f].vertices[0] = i;
				faces[f].vertices[1] = monomers[i].blist[j][0];
				faces[f].vertices[2] = monomers[i].blist[j][1];
				f++;
			}
			if(!check_face(i, monomers[i].blist[j][0], monomers[i].blist[j][2], faces, f)) 
			{
				if(i < Ntype[0]*Nbeads[0])
					faces[f].sphere_id = i/Nbeads[0];
				else
					faces[f].sphere_id = Ntype[0]+(i-Ntype[0]*Nbeads[0])/Nbeads[1];

				faces[f].vertices[0] = i;
				faces[f].vertices[1] = monomers[i].blist[j][0];
				faces[f].vertices[2] = monomers[i].blist[j][2];
				f++;
			}
		}
	}
	printf("%d faces added. %d faces / sphere0  %d faces / sphere1, total faces %d\n", 
			f, Nfaces[0], Nfaces[1], Nfaces[0]*Ntype[0]+Nfaces[1]*Ntype[1]);
	if(f != Nfaces[0]*Ntype[0]+Nfaces[1]*Ntype[1])
		fatal_err("Number of faces does not match", -1);
*/

  //WriteBlist(sphere_pm, monomers, work_dir);  //  201607
  AssignBlist(sphere_pm, monomers, work_dir);
  SetFace(sphere_pm, monomers, faces);

	// Initialize sphere properties 
	for(i=0; i < Nsphere; i++) {
		for(j=0; j < DIMS; j++) {
			spheres[i].com[j]=0.0;
			spheres[i].com0[j] = spheres[i].com[j];
			spheres[i].comold[j] = spheres[i].com[j];
		}
	}

	sphere_props(sphere_pm, spheres, monomers, faces);

	for(i=0; i < Nsphere; i++) 
	{
		spheres[i].disp2 = 0.0;
		spheres[i].dr2 = 0.0;
		for(j=0; j < DIMS; j++) {
			spheres[i].com0[j] = spheres[i].com[j];
			spheres[i].comold[j] = spheres[i].com[j];
		}
	}

	for(type=0; type < NTYPES; type++)  
		for(i=0; i < maxNbeads; i++)
	    free(v[type][i]);
	free(v[0]);
	free(v[1]);

	return 1;
}

// Check if monomer m has been added to the blist of monomer k
int check_bond(int k, int m, struct monomer *monomers)
{
	int i;
	for(i=1; i <= monomers[k].blist[0][0] ; i++) // blist[0][0]: # of bonded monomers
		if(monomers[k].blist[i][0] == m) // blist[i][0], i>0 : the label of a bonded monomers
      return 1;
	return 0;
}

/* check if a face already exists */
int check_face(int a, int b, int c, struct face *faces, int nface)
{
	int i,j,d,check;
	int aa, bb, cc;
	for(i=0; i < nface; i++) 
  {
		aa = faces[i].vertices[0];
		bb = faces[i].vertices[1];
		cc = faces[i].vertices[2];

		if(aa == a)
			if((bb == b && cc == c) || ((bb == c) && (cc == b)))
				return 1;
		if(aa == b)
			if((bb == a && cc == c) || ((bb == c) && (cc == a)))
				return 1;
		if(aa == c)
			if((bb == a && cc == b) || ((bb == b) && (cc == a)))
				return 1;
	}
	return 0;
}

/* return the midpoint between a and b */
int get_midpoint(int a, int b, struct monomer *monomers, int **midpoints)
{
	int i;
	for(i=1; i <= monomers[a].blist[0][0]; i++)
		if(monomers[a].blist[i][0] == b)  /* if a is bonded to b */
			return midpoints[a][i];
	for(i=1; i <= monomers[b].blist[0][0]; i++)
		if(monomers[b].blist[i][0] == a)  /* if b is bonded to a */
			return midpoints[b][i];

	fprintf(stderr, "get_midpoint: monomer %d is not bonded to monomer %d!\n", a, b);
	exit(1);
	return -1;
}

// Determine vertices position and bonding list
void setup(double ***v, int blist[NVERTICES][NBONDS+1][3], double radius[NTYPES])
{
  // Set up positions of vertices
	for(int i=0; i < NTYPES; i++) {
	  double edgelength = 4. * radius[i] / sqrt( 10. + 2.*sqrt(5.) );
		double a = edgelength / 2.;
		double b = edgelength * ( 1.+sqrt(5.) ) / 4.;
		v[i][0][0] = 0.;    v[i][0][1] =  a;    v[i][0][2] =  b;
		v[i][1][0] = 0.;    v[i][1][1] =  a;    v[i][1][2] = -b;
		v[i][2][0] = 0.;    v[i][2][1] = -a;    v[i][2][2] =  b;
		v[i][3][0] = 0.;    v[i][3][1] = -a;    v[i][3][2] = -b;
		v[i][4][0] =  a;    v[i][4][1] =  b;    v[i][4][2] = 0.;
		v[i][5][0] =  a;    v[i][5][1] = -b;    v[i][5][2] = 0.;
		v[i][6][0] = -a;    v[i][6][1] =  b;    v[i][6][2] = 0.;
		v[i][7][0] = -a;    v[i][7][1] = -b;    v[i][7][2] = 0.;
		v[i][8][0] =  b;    v[i][8][1] = 0.;    v[i][8][2] =  a;
		v[i][9][0] = -b;    v[i][9][1] = 0.;    v[i][9][2] =  a;
		v[i][10][0] =  b;   v[i][10][1] = 0.;   v[i][10][2] = -a;
		v[i][11][0] = -b;   v[i][11][1] = 0.;   v[i][11][2] = -a;
	}
	// Set up the bonding list
	/* blist[i][j][0] is the j-th neighbor of vertex i */
	/* blist[i][j][1] is the neighbor of vertex i that makes a 
		 triangle together with blist[i][j][0] */
	/* blist[i][j][2] is another such neighbor (see below) */
	/*                i
  									/|\
	  							 / | \
   blist[i][j][1] /  |  \
	 								\  |  / blist[i][j][2]
									 \ | /
										\|/
				 blist[i][j][0]     */
	/* it is important that with right hand rule, 
		 i-blist[i][j][1]-blist[i][j][0] and i-blist[i][j][0]-blist[i][j][2] 
		 both point toward outside of the sphere!! */
	blist[0][1][0] = 2;   blist[0][1][1] = 9;   blist[0][1][2] = 8;
	blist[0][2][0] = 4;   blist[0][2][1] = 8;   blist[0][2][2] = 6;
	blist[0][3][0] = 6;   blist[0][3][1] = 4;   blist[0][3][2] = 9;
	blist[0][4][0] = 8;   blist[0][4][1] = 2;   blist[0][4][2] = 4;
	blist[0][5][0] = 9;   blist[0][5][1] = 6;   blist[0][5][2] = 2;

	blist[1][1][0] = 3;   blist[1][1][1] = 10;  blist[1][1][2] = 11;
	blist[1][2][0] = 4;   blist[1][2][1] = 6;   blist[1][2][2] = 10;
	blist[1][3][0] = 6;   blist[1][3][1] = 11;  blist[1][3][2] = 4;
	blist[1][4][0] = 10;  blist[1][4][1] = 4;   blist[1][4][2] = 3;
	blist[1][5][0] = 11;  blist[1][5][1] = 3;   blist[1][5][2] = 6;

	blist[2][1][0] = 0;   blist[2][1][1] = 8;   blist[2][1][2] = 9;
	blist[2][2][0] = 5;   blist[2][2][1] = 7;   blist[2][2][2] = 8;
	blist[2][3][0] = 7;   blist[2][3][1] = 9;   blist[2][3][2] = 5;
	blist[2][4][0] = 8;   blist[2][4][1] = 5;   blist[2][4][2] = 0;
	blist[2][5][0] = 9;   blist[2][5][1] = 0;   blist[2][5][2] = 7;

	blist[3][1][0] = 1;   blist[3][1][1] = 11;  blist[3][1][2] = 10;
	blist[3][2][0] = 5;   blist[3][2][1] = 10;  blist[3][2][2] = 7;
	blist[3][3][0] = 7;   blist[3][3][1] = 5;   blist[3][3][2] = 11;
	blist[3][4][0] = 10;  blist[3][4][1] = 1;   blist[3][4][2] = 5;
	blist[3][5][0] = 11;  blist[3][5][1] = 7;   blist[3][5][2] = 1;

	blist[4][1][0] = 0;   blist[4][1][1] = 6;   blist[4][1][2] = 8;
	blist[4][2][0] = 1;   blist[4][2][1] = 10;  blist[4][2][2] = 6;
	blist[4][3][0] = 6;   blist[4][3][1] = 1;   blist[4][3][2] = 0;
	blist[4][4][0] = 8;   blist[4][4][1] = 0;   blist[4][4][2] = 10;
	blist[4][5][0] = 10;  blist[4][5][1] = 8;   blist[4][5][2] = 1;

	blist[5][1][0] = 2;   blist[5][1][1] = 8;   blist[5][1][2] = 7;
	blist[5][2][0] = 3;   blist[5][2][1] = 7;   blist[5][2][2] = 10;
	blist[5][3][0] = 7;   blist[5][3][1] = 2;   blist[5][3][2] = 3;
	blist[5][4][0] = 8;   blist[5][4][1] = 10;  blist[5][4][2] = 2;
	blist[5][5][0] = 10;  blist[5][5][1] = 3;   blist[5][5][2] = 8;

	blist[6][1][0] = 0;   blist[6][1][1] = 9;   blist[6][1][2] = 4;
	blist[6][2][0] = 1;   blist[6][2][1] = 4;   blist[6][2][2] = 11;
	blist[6][3][0] = 4;   blist[6][3][1] = 0;   blist[6][3][2] = 1;
	blist[6][4][0] = 9;   blist[6][4][1] = 11;  blist[6][4][2] = 0;
	blist[6][5][0] = 11;  blist[6][5][1] = 1;   blist[6][5][2] = 9;

	blist[7][1][0] = 2;   blist[7][1][1] = 5;   blist[7][1][2] = 9;
	blist[7][2][0] = 3;   blist[7][2][1] = 11;  blist[7][2][2] = 5;
	blist[7][3][0] = 5;   blist[7][3][1] = 3;   blist[7][3][2] = 2;
	blist[7][4][0] = 9;   blist[7][4][1] = 2;   blist[7][4][2] = 11;
	blist[7][5][0] = 11;  blist[7][5][1] = 9;   blist[7][5][2] = 3;

	blist[8][1][0] = 0;   blist[8][1][1] = 4;   blist[8][1][2] = 2;
	blist[8][2][0] = 2;   blist[8][2][1] = 0;   blist[8][2][2] = 5;
	blist[8][3][0] = 4;   blist[8][3][1] = 10;  blist[8][3][2] = 0;
	blist[8][4][0] = 5;   blist[8][4][1] = 2;   blist[8][4][2] = 10;
	blist[8][5][0] = 10;  blist[8][5][1] = 5;   blist[8][5][2] = 4;

	blist[9][1][0] = 0;   blist[9][1][1] = 2;   blist[9][1][2] = 6;
	blist[9][2][0] = 2;   blist[9][2][1] = 7;   blist[9][2][2] = 0;
	blist[9][3][0] = 6;   blist[9][3][1] = 0;   blist[9][3][2] = 11;
	blist[9][4][0] = 7;   blist[9][4][1] = 11;  blist[9][4][2] = 2;
	blist[9][5][0] = 11;  blist[9][5][1] = 6;   blist[9][5][2] = 7;

	blist[10][1][0] = 1;  blist[10][1][1] = 3;  blist[10][1][2] = 4;
	blist[10][2][0] = 3;  blist[10][2][1] = 5;  blist[10][2][2] = 1;
	blist[10][3][0] = 4;  blist[10][3][1] = 1;  blist[10][3][2] = 8;
	blist[10][4][0] = 5;  blist[10][4][1] = 8;  blist[10][4][2] = 3;
	blist[10][5][0] = 8;  blist[10][5][1] = 4;  blist[10][5][2] = 5;

	blist[11][1][0] = 1;  blist[11][1][1] = 6;  blist[11][1][2] = 3;  
	blist[11][2][0] = 3;  blist[11][2][1] = 1;  blist[11][2][2] = 7;
	blist[11][3][0] = 6;  blist[11][3][1] = 9;  blist[11][3][2] = 1;
	blist[11][4][0] = 7;  blist[11][4][1] = 3;  blist[11][4][2] = 9;
	blist[11][5][0] = 9;  blist[11][5][1] = 7;  blist[11][5][2] = 6;
}

/* determine vertices position and bonding list */
void setup_template(double **v, int ***blist, int nvertices, double radius, char *work_dir)
{
	int i,j;
	double a, b, edgelength;
	char filename[200];
	FILE *stream;

	/* set the coordinates of an icosahedron of the given radius */
//	sprintf(filename, "%s/init/template.init", work_dir);
	sprintf(filename, "%s/init/n2_biconcave.dat", work_dir);
  stream = fopen(filename, "r");

	for(i=0; i < nvertices; i++) {
	  fscanf(stream, "%le %le %le\n", &v[i][0], &v[i][1], &v[i][2]);
    // v[i][0] *= (radius / 4.43);
    // v[i][1] *= (radius / 4.43);
    // v[i][2] *= (radius / 4.43);
		v[i][0] *= 0.5;
		v[i][1] *= 0.5; 	
		v[i][2] *= 0.5;
  }
	for(i=0; i < nvertices; i++) {
		fscanf(stream, "%*s %d\n", &blist[i][0][0]);
		for(j=1; j<=blist[i][0][0]; j++)
			fscanf(stream, "%d %d %d", &blist[i][j][0], &blist[i][j][1], &blist[i][j][2]);
	}
	fclose(stream);
}

void write_particle_para(struct sphere_param *sphere_pm, struct monomer *vertex, struct 
     face *faces, char *work_dir)
{
  int numVertex = sphere_pm->N_per_sphere[0];
  for(int n1=0; n1 < numVertex; n1++)
  {
    for(int label=1; label <= vertex[n1].blist[0][0]; label++)
    {
      // calculate initial bond length
      int n2 = vertex[n1].blist[label][0];
      double bond[DIMS];
      double bondLength2 = 0.;
      for(int d=0; d < DIMS; d++) {
        bond[d] = vertex[n1].pos_pbc[d] - vertex[n2].pos_pbc[d];
        bondLength2 += bond[d]*bond[d];
      }
      vertex[n1].initLength[label] = sqrt(bondLength2);

      // calculate initial bending angle 
      int n3 = vertex[n1].blist[label][1];
      int n4 = vertex[n1].blist[label][2];
      double x31[DIMS], x21[DIMS], x41[DIMS];
      double normal1[DIMS], normal2[DIMS];
      double cross[DIMS], crossNorm;
      for(int d=0; d < DIMS; d++) {
        x31[d] = vertex[n3].pos_pbc[d] - vertex[n1].pos_pbc[d];
        x21[d] = vertex[n2].pos_pbc[d] - vertex[n1].pos_pbc[d];
        x41[d] = vertex[n4].pos_pbc[d] - vertex[n1].pos_pbc[d];
      }
      product(x31,x21, normal1);
      product(x21,x41, normal2);
      product(normal1,normal2, cross);
      crossNorm = cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2];
      crossNorm = sqrt(crossNorm);
      vertex[n1].initAngle[label] = atan2(crossNorm,iproduct(normal1,normal2));
    }
  }
  // calculate the center of mass
  double com[3];
  for(int d=0; d < DIMS; d++) {
    com[d] = 0.0;
    for(int j=0; j < numVertex; j++)
      com[d] += vertex[j].pos_pbc[d];
    com[d] /= numVertex;
  }
  // calculate initial area and volume
  double area=0.;
  double volume=0.;
  for(int nface=0; nface < sphere_pm->face_per_sphere[0]; nface++)
  {
    int n1 = faces[nface].vertices[0];
    int n2 = faces[nface].vertices[1];
    int n3 = faces[nface].vertices[2];
    double dr[3], q1[3], q2[3], normal[3];
    for(int d=0; d < DIMS; d++)
      dr[d] = vertex[n1].pos_pbc[d] - com[d];
    for(int d=0; d < DIMS; d++)
      q1[d] = vertex[n2].pos_pbc[d] - vertex[n1].pos_pbc[d];
    for(int d=0; d < DIMS; d++)
      q2[d] = vertex[n3].pos_pbc[d] - vertex[n1].pos_pbc[d];
    product(q1, q2, normal);
    faces[nface].area_0 = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2])/2.0;
    volume += fabs(dr[0]*normal[0]+dr[1]*normal[1]+dr[2]*normal[2])/6.;
    area += faces[nface].area_0;
  }
  sphere_pm->V0[0] = volume;
  sphere_pm->A0[0] = area;
  //printf("V0=%le  A0=%le\n", volume, area);
  // Write data on the disk
  char fileName[100];
  FILE *stream;
  sprintf(fileName, "%s/init/initShapePara_n2.dat",work_dir);
  stream = fopen(fileName,"w");
  fprintf(stream,"%le %le\n",volume, area);
  for(int n1=0; n1 < numVertex; n1++)
    for(int label=1; label <= vertex[n1].blist[0][0]; label++)
      fprintf(stream,"%le  %le\n",vertex[n1].initLength[label], 
              vertex[n1].initAngle[label]);
  for(int nface=0; nface < sphere_pm->face_per_sphere[0]; nface++)
    fprintf(stream,"%le\n",faces[nface].area_0);
  fclose(stream);
}

void calculate_particle_para(struct sphere_param *sphere_pm, struct monomer *vertex, 
     struct face *faces, char *work_dir)
{
  int numVertex = sphere_pm->N_per_sphere[0];
  for(int n1=0; n1 < numVertex; n1++)
  {
    for(int label=1; label <= vertex[n1].blist[0][0]; label++)
    {
      // calculate initial bond length
      int n2 = vertex[n1].blist[label][0];
      double bond[DIMS];
      double bondLength2 = 0.;
      for(int d=0; d < DIMS; d++) {
        bond[d] = vertex[n1].pos_pbc[d] - vertex[n2].pos_pbc[d];
        bondLength2 += bond[d]*bond[d];
      }
      vertex[n1].initLength_temp[label] = sqrt(bondLength2); 

      // calculate initial bending angle 
      int n3 = vertex[n1].blist[label][1];
      int n4 = vertex[n1].blist[label][2];
      double x31[DIMS], x21[DIMS], x41[DIMS];
      double normal1[DIMS], normal2[DIMS];
      double cross[DIMS], crossNorm;
      for(int d=0; d < DIMS; d++) {
        x31[d] = vertex[n3].pos_pbc[d] - vertex[n1].pos_pbc[d];
        x21[d] = vertex[n2].pos_pbc[d] - vertex[n1].pos_pbc[d];
        x41[d] = vertex[n4].pos_pbc[d] - vertex[n1].pos_pbc[d];
      }
      product(x31,x21, normal1);
      product(x21,x41, normal2);
      product(normal1,normal2, cross);
      crossNorm = cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2];
      crossNorm = sqrt(crossNorm);
      vertex[n1].initAngle_temp[label] = atan2(crossNorm,iproduct(normal1,normal2));
    }
  }
  // calculate the center of mass
  double com[3];
  for(int d=0; d < DIMS; d++) {
    com[d] = 0.0;
    for(int j=0; j < numVertex; j++)
      com[d] += vertex[j].pos_pbc[d];
    com[d] /= numVertex;
  }
  // calculate initial area and volume
  double area=0.;
  double volume=0.;
  for(int nface=0; nface < sphere_pm->face_per_sphere[0]; nface++)
  {
    int n1 = faces[nface].vertices[0];
    int n2 = faces[nface].vertices[1];
    int n3 = faces[nface].vertices[2];
    double dr[3], q1[3], q2[3], normal[3];
    for(int d=0; d < DIMS; d++)
      dr[d] = vertex[n1].pos_pbc[d] - com[d];
    for(int d=0; d < DIMS; d++)
      q1[d] = vertex[n2].pos_pbc[d] - vertex[n1].pos_pbc[d];
    for(int d=0; d < DIMS; d++)
      q2[d] = vertex[n3].pos_pbc[d] - vertex[n1].pos_pbc[d];
    product(q1, q2, normal);
    volume += fabs(dr[0]*normal[0]+dr[1]*normal[1]+dr[2]*normal[2])/6.;
    area += sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2])/2.0;
  }
  sphere_pm->V0_temp[0] = volume;
  sphere_pm->A0_temp[0] = area;

  // write data on the disk
  printf("At line 1714 in init_sphere.c, it's unnecessary to write data on the disk!\n");
  char fileName[100];
  FILE *stream;
  sprintf(fileName, "%s/init/particlePara_temp_n2.dat",work_dir);
  stream = fopen(fileName,"w");
  fprintf(stream,"%le %le\n",volume, area);
  for(int n1=0; n1 < numVertex; n1++)
    for(int label=1; label <= vertex[n1].blist[0][0]; label++)
      fprintf(stream,"%le  %le\n",vertex[n1].initLength_temp[label], 
              vertex[n1].initAngle_temp[label]);
  fclose(stream);
}

void read_particle_para(struct sphere_param *sphere_pm, struct monomer *vertex, struct 
     face *faces, char *work_dir)
{
  char fileName[150];
  FILE *stream; 
  int numVertex = sphere_pm->N_per_sphere[0];
  int numParticle = sphere_pm->Ntype[0];
  sprintf(fileName,"%s/init/initShapePara_n2.dat",work_dir);
  stream = fopen(fileName,"r");

  // For the 1st particle
  fscanf(stream,"%le %le\n", &sphere_pm->V0[0], &sphere_pm->A0[0]);
  for(int n1=0; n1 < numVertex; n1++)
    for(int label=1; label <= vertex[n1].blist[0][0]; label++)
      fscanf(stream, "%le %le\n", &vertex[n1].initLength[label], 
             &vertex[n1].initAngle[label]);
  for(int nface=0; nface < sphere_pm->face_per_sphere[0]; nface++)
    fscanf(stream,"%le\n",&faces[nface].area_0);
  fclose(stream);

  for(int par=1; par < numParticle; par++) // For other particles
  {
    int initVertexLabel = par * numVertex;
    int initFaceLabel = par * sphere_pm->face_per_sphere[0];

    for(int j=0; j < numVertex; j++)
    {
      int n1 = j + initVertexLabel;
      for(int bond=1; bond <= vertex[n1].blist[0][0]; bond++) {
        vertex[n1].initLength[bond] = vertex[j].initLength[bond];
        vertex[n1].initAngle[bond] = vertex[j].initAngle[bond];
//printf("vertex[%d].initLength[%d]=%f\nvertex[%d].initLength[%d]=%f\n",
//n1,bond,vertex[n1].initLength[bond],j,bond,vertex[j].initLength[bond]);
      }
    }
    for(int nface=0; nface < sphere_pm->face_per_sphere[0]; nface++) 
      faces[initFaceLabel+nface].area_0 = faces[nface].area_0;
  }

  sphere_pm->V0_final[0] = sphere_pm->V0[0];
  sphere_pm->A0_final[0] = sphere_pm->A0[0];
  for(int par=0; par < sphere_pm->Ntype[0]; par++) {
    int initVertexLabel = par * sphere_pm->N_per_sphere[0];
    for(int j=0; j < sphere_pm->N_per_sphere[0]; j++) {
      int n1 = j + initVertexLabel;
      for(int bond=1; bond <= vertex[n1].blist[0][0]; bond++) {
        vertex[n1].initLength_final[bond] = vertex[n1].initLength[bond];
      }
    }
  }
}

void biconcaveTemplate(struct sphere_param *sphere_pm, struct monomer *vertex)
{
  // calculate the center of mass
  double com[3];
  int numVertex = sphere_pm->N_per_sphere[0];
  for(int d=0; d < DIMS; d++) {
    com[d] = 0.0;
    for(int j=0; j < numVertex; j++)
      com[d] += vertex[j].pos_pbc[d];
    com[d] /= numVertex;
  }
  for(int n1=0; n1 < numVertex; n1++) {
    vertex[n1].pos_com[0] = vertex[n1].pos_pbc[0] - com[0];    
    vertex[n1].pos_com[1] = vertex[n1].pos_pbc[1] - com[1];    
    vertex[n1].pos_com[2] = vertex[n1].pos_pbc[2] - com[2];
  } 

  double c0 = 0.0518;
  double c2 = 2.0026;
  double c4 = -4.491;
  double r = 3.91;
	extern int max_x, max_y, max_z;
  double maxsize[3];
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;
  for(int n1=0; n1 < numVertex; n1++)
  {
    double rho = vertex[n1].pos_com[0] * vertex[n1].pos_com[0] + vertex[n1].pos_com[2] *
                 vertex[n1].pos_com[2];
    rho = sqrt(rho);
    double rho_over_r = rho / r;
    double rho_over_r2 = rho_over_r * rho_over_r;
    double rho_over_r4 = rho_over_r2 * rho_over_r2;
    if(vertex[n1].pos_com[1] >= 0)
      vertex[n1].pos_com[1] = 2*r * sqrt(1 - rho_over_r2) * (c0 + c2/4*rho_over_r2 + 
                                   c4/16*rho_over_r4);
    else
      vertex[n1].pos_com[1] = -2*r * sqrt(1 - rho_over_r2) * (c0 + c2/4*rho_over_r2 + 
                                    c4/16*rho_over_r4);

    vertex[n1].pos_pbc[0] = vertex[n1].pos_com[0] + com[0];
    vertex[n1].pos_pbc[1] = vertex[n1].pos_com[1] + com[1];
    vertex[n1].pos_pbc[2] = vertex[n1].pos_com[2] + com[2];
    vertex[n1].pos[0] = box(vertex[n1].pos_pbc[0], maxsize[0]);
    vertex[n1].pos[1] = box(vertex[n1].pos_pbc[1], maxsize[1]);
    vertex[n1].pos[2] = box(vertex[n1].pos_pbc[2], maxsize[2]);
  }
}

void growth_procedure(int n_step, int totalGrowthStep, struct sphere_param *sphere_pm, 
            struct monomer *monomers, struct object *objects, int n_cycle, int num_step,
            int mark_interval, int oscillation_period, int window, struct sphere *spheres,
            struct face *faces, Float ***velcs_df, int **node_map, char *work_dir, 
            VSLStreamStatePtr rngstream)
{
  extern int num_sph;
  extern int    max_x, max_y, max_z;
  static struct vector f_tot, p_lbe;
  double sum_y[MAX_Y][16];
  FILE *file_ptr=0;

  int type, numVertex, initVertexLabel, numParticle, Nsphere[NTYPES];
  Nsphere[0]=sphere_pm->Ntype[0];
  Nsphere[1]=sphere_pm->Ntype[1];
  numParticle = Nsphere[0] + Nsphere[1];

  if(n_step < totalGrowthStep) 
  {
    sphere_pm->V0[0] = sphere_pm->V0_temp[0] + 
            n_step * (sphere_pm->V0_final[0]- sphere_pm->V0_temp[0]) / totalGrowthStep;
    sphere_pm->A0[0] = sphere_pm->A0_temp[0] + 
            n_step * (sphere_pm->A0_final[0]- sphere_pm->A0_temp[0]) / totalGrowthStep;
    for(int par=0; par < numParticle; par++)
    {
      type = (par < sphere_pm->Ntype[0] ? 0:1);
 	    numVertex = sphere_pm->N_per_sphere[type];
 		  initVertexLabel = (type==0 ? par*sphere_pm->N_per_sphere[0]: 
      Nsphere[0]*sphere_pm->N_per_sphere[0] + (par-Nsphere[0])*sphere_pm->N_per_sphere[1]);

      for(int j=0; j < numVertex; j++) {
        int n1 = initVertexLabel + j;
        for(int bond=1; bond <= monomers[n1].blist[0][0]; bond++) 
          monomers[n1].initLength[bond] = monomers[j].initLength_temp[bond] + 
          n_step * (monomers[j].initLength_final[bond]-monomers[j].initLength_temp[bond]) /
          totalGrowthStep;        
      }
    }
  }
  else 
  {
    //sphere_pm->k_V[0]=10000.0;
    //sphere_pm->k_A[0]=20000.0; // modification 20170314 temporary design! 
    sphere_pm->V0[0] = sphere_pm->V0_final[0];
    sphere_pm->A0[0] = sphere_pm->A0_final[0];
    for(int par=0; par < numParticle; par++)
    {
      type = (par < sphere_pm->Ntype[0] ? 0:1);
 		  numVertex = sphere_pm->N_per_sphere[type];
 		  initVertexLabel = (type==0 ? par*sphere_pm->N_per_sphere[0]: 
                        Nsphere[0]*sphere_pm->N_per_sphere[0] + 
                        (par-Nsphere[0])*sphere_pm->N_per_sphere[1]);
      for(int j=0; j < numVertex; j++) {
        int n1 = initVertexLabel + j;
        for(int bond=1; bond <= monomers[n1].blist[0][0]; bond++) 
          monomers[n1].initLength[bond] = monomers[j].initLength_final[bond];    
      }  
    }
  } 
  if(n_step==10000) {
    sphere_pm->k_V[0]=10000.0;
    sphere_pm->k_A[0]=20000.0;
    objects[0].u = objects[0].u_init;
    objects[1].u = objects[1].u_init;
  }

  double dr2_temp = 0.0;
  double dr2[2] = {0.0, 0.0};
  double dx, dy, dz;
	for(int n_sph = 0; n_sph < sphere_pm->num_beads; n_sph++)
	{
	  dx = n_image(monomers[n_sph].pos[0] - monomers[n_sph].pos_lst[0], max_x);
		dy = n_image(monomers[n_sph].pos[1] - monomers[n_sph].pos_lst[1], max_y);
		dz = n_image(monomers[n_sph].pos[2] - monomers[n_sph].pos_lst[2], max_z);
    dr2_temp = dx*dx + dy*dy + dz*dz;

    if(dr2_temp > dr2[1] && dr2_temp > dr2[0]) {
      dr2[1] = dr2[0];
      dr2[0] = dr2_temp;
    }
    else if(dr2_temp >= dr2[1] && dr2_temp < dr2[0]) 
      dr2[1] = dr2_temp;
  }
  if((sqrt(dr2[0]) + sqrt(dr2[1])) > 1e-3) { // Modification 20170316 0.5->0.25
	  n_list_mon(sphere_pm, monomers);  // Update neighbor lists
	}

  Write_Output (n_step, n_cycle, num_step, mark_interval, oscillation_period, window, 
                num_sph, objects, spheres, monomers, faces, sphere_pm,velcs_df, 
                node_map, work_dir);

  double fluidStress_pre[6]={0,0};
  double fluidStress_pos[6]={0,0};
  lbe_update (objects, sphere_pm, spheres, monomers, velcs_df, node_map, f_tot, &p_lbe,
      sum_y, num_step, n_step, file_ptr, rngstream, fluidStress_pre,
      fluidStress_pos);

  for(int step=0; step < sphere_pm->MD_steps; step++)
    verlet_update (monomers, faces, sphere_pm, velcs_df, n_step, rngstream); 
}

void growth_procedure_2(int n_step, int totalGrowthStep, struct sphere_param *sphere_pm, 
            struct monomer *monomers, struct object *objects, int n_cycle, int num_step,
            int mark_interval, int oscillation_period, int window, struct sphere *spheres,
            struct face *faces, Float ***velcs_df, int **node_map, char *work_dir, 
            VSLStreamStatePtr rngstream)
{
  extern int num_sph, max_x, max_y, max_z;
  int type, numVertex, initVertexLabel, numParticle, Nsphere[NTYPES];
  Nsphere[0]=sphere_pm->Ntype[0];
  Nsphere[1]=sphere_pm->Ntype[1];
  numParticle = Nsphere[0] + Nsphere[1];

  // Restore V0, A0, and spring length gradually
  if(n_step < totalGrowthStep) 
  {
    sphere_pm->V0[0] = sphere_pm->V0_temp[0] + 
    n_step * (sphere_pm->V0_final[0]- sphere_pm->V0_temp[0]) / totalGrowthStep;
    sphere_pm->A0[0] = sphere_pm->A0_temp[0] + 
    n_step * (sphere_pm->A0_final[0]- sphere_pm->A0_temp[0]) / totalGrowthStep;
    for(int par=0; par < numParticle; par++)
    {
      type = (par < sphere_pm->Ntype[0] ? 0:1);
 	    numVertex = sphere_pm->N_per_sphere[type];
 		  initVertexLabel = (type==0 ? par*sphere_pm->N_per_sphere[0]: 
      Nsphere[0]*sphere_pm->N_per_sphere[0] + (par-Nsphere[0])*sphere_pm->N_per_sphere[1]);

      for(int j=0; j < numVertex; j++) {
        int n1 = initVertexLabel + j;
        for(int bond=1; bond <= monomers[n1].blist[0][0]; bond++) 
          monomers[n1].initLength[bond] = monomers[j].initLength_temp[bond] + 
          n_step * (monomers[j].initLength_final[bond]-monomers[j].initLength_temp[bond]) /
          totalGrowthStep;        
      }
    }

    sphere_pm->k_V[0] =   (double) (n_step / totalGrowthStep) * 50.0;
    sphere_pm->k_A[0] =   (double) (n_step / totalGrowthStep) * 100.0;
    sphere_pm->ka_local = (double) (n_step / totalGrowthStep) * 10.0;

  }
  else 
  {
    sphere_pm->V0[0] = sphere_pm->V0_final[0];
    sphere_pm->A0[0] = sphere_pm->A0_final[0];
    for(int par=0; par < numParticle; par++)
    {
      type = (par < sphere_pm->Ntype[0] ? 0:1);
 		  numVertex = sphere_pm->N_per_sphere[type];
 		  initVertexLabel = (type==0 ? par*sphere_pm->N_per_sphere[0]: 
                        Nsphere[0]*sphere_pm->N_per_sphere[0] + 
                        (par-Nsphere[0])*sphere_pm->N_per_sphere[1]);
      for(int j=0; j < numVertex; j++) {
        int n1 = initVertexLabel + j;
        for(int bond=1; bond <= monomers[n1].blist[0][0]; bond++) 
          monomers[n1].initLength[bond] = monomers[j].initLength_final[bond];    
      }  
    }
    sphere_pm->k_V[0] = 50.0;
    sphere_pm->k_A[0] = 100.0;
    sphere_pm->ka_local = 10.0;
  } 
  //// Apply area and volume restoring forces gradually
  //extern int initStep; // Modification 20170712: Temporary design
  //double step_temp = totalGrowthStep + initStep;
  //if(n_step > totalGrowthStep && n_step <= step_temp) {
  //  sphere_pm->k_V[0] = (n_step-totalGrowthStep) / (step_temp-totalGrowthStep) * 50.0;
  //  sphere_pm->k_A[0] = (n_step-totalGrowthStep) / (step_temp-totalGrowthStep) * 100.0;
  //  sphere_pm->ka_local = (n_step-totalGrowthStep) / (step_temp-totalGrowthStep) * 10.0;
  //}
  // Turn on the simple shear flow
  if(n_step == totalGrowthStep) {
    objects[0].u = objects[0].u_init;
    objects[1].u = objects[1].u_init;
//    sphere_pm->attrac_type = sphere_pm->attrac_type_1;  // Modification 20170802: Temporary!!!
  }
  // Turn on adhesive interaction
  // Modification 20170711
  // Temporary design: Turn on local area constraint force
  //if(n_step == step_temp)
  //{
  //  sphere_pm->attrac_type = 2;
  //  //sphere_pm->ka_local = 10.;
  //}
  // Update the neighbor list 
  double dr2_temp = 0.0;
  double dr2[2] = {0.0, 0.0};
  double dx, dy, dz;
	for(int n_sph = 0; n_sph < sphere_pm->num_beads; n_sph++)
	{
	  dx = n_image(monomers[n_sph].pos[0] - monomers[n_sph].pos_lst[0], max_x);
		dy = n_image(monomers[n_sph].pos[1] - monomers[n_sph].pos_lst[1], max_y);
		dz = n_image(monomers[n_sph].pos[2] - monomers[n_sph].pos_lst[2], max_z);
    dr2_temp = dx*dx + dy*dy + dz*dz;

    if(dr2_temp > dr2[1] && dr2_temp > dr2[0]) {
      dr2[1] = dr2[0];
      dr2[0] = dr2_temp;
    }
    else if(dr2_temp >= dr2[1] && dr2_temp < dr2[0]) 
      dr2[1] = dr2_temp;
  }
  if((sqrt(dr2[0]) + sqrt(dr2[1])) > 1e-3) { // Modification 20170316 0.5->0.25
	  n_list_mon(sphere_pm, monomers);  // Update neighbor lists
	}

  //Write_Output (n_step, n_cycle, num_step, mark_interval, oscillation_period, window, 
  //    num_sph,objects, spheres, monomers, faces, sphere_pm,velcs_df, node_map, work_dir);
	if(n_step % 1000 == 0)  
    Write_time (n_step,work_dir); 
  Write_Monomer (spheres,monomers,sphere_pm,n_step,n_cycle*num_step,work_dir); 
	Write_Sphere (spheres,monomers,faces,sphere_pm,n_step,n_cycle*num_step,work_dir);
	write_growth_config (spheres,monomers,faces,sphere_pm,n_step,n_cycle*num_step,work_dir);
	if(n_step % sphere_pm->write_fluid == 0 && sphere_pm->write_fluid < num_step) 
		Write_Fluid (velcs_df,node_map,n_cycle*num_step+n_step,work_dir);

  get_force_growth (sphere_pm, monomers, faces, velcs_df, n_step, rngstream);
  euler_update (sphere_pm->num_beads, monomers);
}

