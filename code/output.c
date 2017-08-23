#include "header.h"
double springLength_eq=0.;  // modified 20160914
void sym33eigen_values(double mx[3][3], double lambda[3]);
void eigen_vector(double mx[3][3], double lambda, double v[3]);

void Write_Output (int n_step, int n_cycle, int num_step, int mark_interval, 
     int oscillation_period, int window, int num_sph, struct object *objects, 
     struct sphere *spheres, struct monomer *monomers, struct face *faces, struct 
     sphere_param *sphere_pm, Float ***velcs_df, int **node_map, char *work_dir)
{
	if(n_step % 1000 == 0)  
    Write_time (n_step+n_cycle*num_step, work_dir); 
	if(sphere_pm->num_beads !=0) {
		Write_Monomer (spheres, monomers, sphere_pm, n_step, n_cycle*num_step, work_dir);  // Modification 20170316
		Write_Sphere (spheres, monomers, faces, sphere_pm, n_step, n_cycle*num_step, work_dir);
		Write_Sphereconfig (spheres, monomers, faces, sphere_pm, n_step, n_cycle*num_step, 
				                work_dir);
	}

	if(n_step % sphere_pm->write_fluid == 0 && sphere_pm->write_fluid < num_step) {
		Write_Fluid (velcs_df, node_map, n_cycle*num_step+n_step, work_dir);
		//    Write_Nodemap(node_map,  n_cycle*num_step+n_step, work_dir);
//		Write_Velslice (monomers, velcs_df, node_map, n_cycle*num_step+n_step, work_dir); 
	}
}

void sphere_props(struct sphere_param *sphere_pm, struct sphere *spheres, struct monomer 
		 *mono, struct face *faces)
{
	extern int max_x, max_y, max_z;
	int i, j, k, d, n1, n2, n3, d1, d2;
	int Nsphere = sphere_pm->Nsphere;
	int Ntype[NTYPES], Nbead[NTYPES];
	int maxsize[3];
	int nspring, nface;
	int start, type;
	int min;
	double sprng_len, avg_sprng_len;
	double temp, tempstretch;
	double g_tensor[9];  /* gyration tensor */
	//double I_tensor[DIMS][DIMS];  /* inertia tensor */
	double dr[DIMS], lambda[DIMS], Im[DIMS];
	double q12[DIMS], q13[DIMS], normal[DIMS];
	double y_com;
	double Rx2, Ry2, Rz2;
	double costheta, sintheta;
	double omega[DIMS];
	double y[162];
	double v[DIMS];  
	FILE *stream;
	double disp2_temp = 0.0;
	double dr2_temp = 0.0;
	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;
	for (i=0 ; i<NTYPES ; i++) {
		Ntype[i] = sphere_pm->Ntype[i];
		Nbead[i] = sphere_pm->N_per_sphere[i];
	}

	// Calculate the sphere center-of-mass
	for(i=0; i < Nsphere; i++) 
	{
		type = (i < Ntype[0] ? 0 : 1);
		start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0] + (i-Ntype[0])*Nbead[1]);

		spheres[i].com[0] = 0.0;
		spheres[i].com[1] = 0.0;
		spheres[i].com[2] = 0.0;
		spheres[i].disp2 = 0.0;  // modified by ctliao 10/14/2015
		spheres[i].dr2 = 0.0;

		for(j=0; j < Nbead[type]; j++) 
			for(d=0; d < DIMS; d++)
				spheres[i].com[d] += mono[start+j].pos_pbc[d];

		for(d=0; d < DIMS; d++) // modified by ctliao 10/14/2015
		{
			spheres[i].com[d] = spheres[i].com[d] / Nbead[type];
			//disp2_temp = spheres[i].com[d] - spheres[i].com0[d];
			//dr2_temp = spheres[i].com[d]-spheres[i].comold[d];
			//spheres[i].disp2 += disp2_temp * disp2_temp;
			//spheres[i].dr2 += dr2_temp * dr2_temp;
			spheres[i].comold[d] = spheres[i].com[d];
		}
	}
	// Calculate the sphere radius of gyration and stretch
	for(i=0; i < Nsphere; i++) 
	{
    double eigenvector[9];
    int it_num = 0;
    int rot_num = 0;
		type = (i < Ntype[0] ? 0 : 1);
		start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0] + (i-Ntype[0])*Nbead[1]);
	  spheres[i].stretch[0] = 0.0;
	  spheres[i].stretch[1] = 0.0;
	  spheres[i].stretch[2] = 0.0;
		for(d=0; d < 9; d++)
		  g_tensor[d] = 0.0;

		for(j=0; j < Nbead[type]; j++) 
		{
      // Calculate the gyration tensor
			for(d=0; d < DIMS; d++)
				dr[d] = mono[start+j].pos_pbc[d] - spheres[i].com[d];
			for(d1=0 ; d1 < DIMS ; d1++)
				for(d2=0 ; d2 < DIMS ; d2++)
					g_tensor[d1*DIMS+d2] += dr[d1]*dr[d2];// {00,01,02,10,11,12,20,21,22}

      // Calculate the degree of stretching
			for(k=0; k < Nbead[type]; k++) 
      {
				for(d=0; d < DIMS; d++) {
					temp = mono[start+k].pos_pbc[d] - mono[start+j].pos_pbc[d];
					if(temp > spheres[i].stretch[d])
						spheres[i].stretch[d] = temp;
				}
			}
		}
		for(d1=0; d1 < 9; d1++)
	    g_tensor[d1] /= Nbead[type];

    // Calculate eigenvalues of the gyration tensor
    jacobi_eigenvalue(3,g_tensor,100,eigenvector,spheres[i].g_eigenvalue,&it_num,
                      &rot_num);
	}
	// Calculate average angular velocity (not sure whether this data useful?)
	for(i=0; i < Nsphere; i++) 
  {
		type = (i < Ntype[0] ? 0 : 1);
		start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0] + (i-Ntype[0])*Nbead[1]);
		for(d=0; d < DIMS; d++)
			spheres[i].omega[d] = 0.;
		for(j=0; j < Nbead[type]; j++) 
    {
			temp = 0.;
			for(d=0; d < DIMS; d++) {
				dr[d] = mono[start+j].pos_pbc[d] - spheres[i].com[d];
				temp += dr[d]*dr[d];
			}
			product(dr, mono[start+j].vel, omega);
			for(d=0; d < DIMS; d++)
				spheres[i].omega[d] += omega[d]/temp;
		}
		for(d=0; d < DIMS; d++)
			spheres[i].omega[d] /= Nbead[type];
	}
	// Calculate the average spring length
	for(i=0; i < Nsphere; i++) 
  {
		type = (i < Ntype[0]? 0 : 1);
		start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0] + (i-Ntype[0])*Nbead[1]);
		nspring = 0;
		avg_sprng_len = 0.0;
		for(j=0; j < Nbead[type]; j++) {
			for(k=1; k <= mono[start+j].blist[0][0]; k++) 
      {
				n1 = start + j;
				n2 = mono[n1].blist[k][0];
				sprng_len = 0.0;
				for(d=0; d < DIMS; d++) {
					temp = n_image(mono[n1].pos[d]-mono[n2].pos[d], maxsize[d]);
					sprng_len += temp*temp;
				}
				sprng_len = sqrt(sprng_len);
				avg_sprng_len += sprng_len;
				nspring ++;
			}
    }
		avg_sprng_len /= nspring;
		spheres[i].avg_sprng_len = avg_sprng_len;
	}
	// Calculate the surface area and volume of each particle
	for(i=0; i < Nsphere; i++) 
  {
		type = (i < Ntype[0] ? 0 : 1);
		start = (type == 0 ? i*sphere_pm->face_per_sphere[0] : 
            Ntype[0]*sphere_pm->face_per_sphere[0] + 
            (i-Ntype[0])*sphere_pm->face_per_sphere[1]);
		spheres[i].area = 0.;
		spheres[i].volume =0.;
		//spheres[i].volume_dupin =0.;
    double inertia[9]={0.};
    double eigenvalue[DIMS]={0.};
    int it_num=0;
    int rot_num=0;
		/* n1, n2, n3 are the three points of the triangle; 
       q12, q13 are two sides; dr is vector from com to n1 */
		for(j=0; j < sphere_pm->face_per_sphere[type]; j++) 
    {
      double unitNormal[DIMS];
      double faceCen[DIMS]; 
      double faceCen_com[DIMS];
      double faceCen_com_mag2=0.;
      
			nface = start + j;
			n1 = faces[nface].vertices[0];
      n2 = faces[nface].vertices[1];
      n3 = faces[nface].vertices[2];
			for(d=0; d < DIMS; d++)
			  dr[d] = mono[n1].pos_pbc[d] - spheres[i].com[d];
			for(d=0; d < DIMS; d++)
				q12[d] = mono[n2].pos_pbc[d] - mono[n1].pos_pbc[d];
			for(d=0; d < DIMS ; d++)
				q13[d] = mono[n3].pos_pbc[d] - mono[n1].pos_pbc[d];
			product(q13, q12, normal);  // "normal" points outward viewing from the centroid of the particle
      double normalMag = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      for(d=0; d < DIMS; d++)
        unitNormal[d] = normal[d] / normalMag;
			faces[nface].area = 0.5 * normalMag;
			spheres[i].area += faces[nface].area;
			spheres[i].volume += fabs(dr[0]*normal[0]+dr[1]*normal[1]+dr[2]*normal[2])/6.;

			//y_com = (mono[n1].pos_pbc[1]+mono[n2].pos_pbc[1]+mono[n3].pos_pbc[1])/3.;
			///* make the y-component of the normal vector point outwards */
			//if(y_com - spheres[i].com[1] >= 0.0) 
			//	spheres[i].volume_dupin += fabs(normal[1]*y_com/2.);  
			//else
			//	spheres[i].volume_dupin -= fabs(normal[1]*y_com/2.); 

      for(d=0; d < DIMS; d++) {
        //faceCen[d] = (mono[n1].pos[d] + mono[n2].pos[d] + mono[n3].pos[d]) / 3.;
        faceCen[d] = (mono[n1].pos_pbc[d] + mono[n2].pos_pbc[d] + mono[n3].pos_pbc[d])/3.; // 20170131
        faceCen_com[d] = faceCen[d] - spheres[i].com[d];
        faceCen_com_mag2 += faceCen_com[d]*faceCen_com[d];
      }
      double factor = iproduct(faceCen_com, unitNormal);
      //0:00, 1:01, 2:02, 3:10, 4:11, 5:12, 6:20, 7:21; 8:22
      inertia[0] += faces[nface].area * (faceCen_com_mag2 - faceCen_com[0]*faceCen_com[0])*factor;
      inertia[1] += faces[nface].area * (-faceCen_com[0]*faceCen_com[1]) * factor;
      inertia[2] += faces[nface].area * (-faceCen_com[0]*faceCen_com[2]) * factor;
      inertia[4] += faces[nface].area * (faceCen_com_mag2 - faceCen_com[1]*faceCen_com[1])*factor;
      inertia[5] += faces[nface].area * (-faceCen_com[1]*faceCen_com[2]) * factor;
      inertia[8] += faces[nface].area * (faceCen_com_mag2 - faceCen_com[2]*faceCen_com[2])*factor;
		}
    inertia[0] *= 0.2; // rho/5 = 0.2, par. density rho divided by 5. 
    inertia[1] *= 0.2;
    inertia[2] *= 0.2;
    inertia[4] *= 0.2;
    inertia[5] *= 0.2;
    inertia[8] *= 0.2;
    inertia[3] = inertia[1];
    inertia[6] = inertia[2];
    inertia[7] = inertia[5];
/*
    double test[16]; 
    double ivector[16]={0.};
    double ivalue[4]={0.};
    test[0]= 4.0;
    test[1]=0.0;
    test[2]=0.0;
    test[3]=0.0;

    test[4]=0.0;
    test[5]=1.0;
    test[6]=0.0;
    test[7]=0.0;

    test[8]=0.0;
    test[9]=0.0;
    test[10]=3.0;
    test[11]=0.0;

    test[12]=0.0;
    test[13]=0.0;
    test[14]=0.0;
    test[15]=2.0;

    jacobi_eigenvalue(4, test, 100, ivector, ivalue, &it_num, &rot_num);
    printf("ivalue= %f %f %f %f\n",ivalue[0],ivalue[1],ivalue[2],ivalue[3]);
    printf("ivector= %f %f %f %f\n",ivector[0],ivector[1],ivector[2],ivector[3]);
    printf("         %f %f %f %f\n",ivector[4],ivector[5],ivector[6],ivector[7]);
    printf("         %f %f %f %f\n",ivector[8],ivector[9],ivector[10],ivector[11]);
    printf("         %f %f %f %f\n",ivector[12],ivector[13],ivector[14],ivector[15]);
*/    
    // Solve the eigenvalue problem of the inertia tensor.  
    jacobi_eigenvalue(3,inertia,100,spheres[i].i_eigenvector,eigenvalue,&it_num,
                      &rot_num);
    spheres[i].semiaxis_a = sqrt(5*(eigenvalue[1] + eigenvalue[2] - eigenvalue[0]) / 
               (2*spheres[i].volume));
    spheres[i].semiaxis_b = sqrt(5*(eigenvalue[2] + eigenvalue[0] - eigenvalue[1]) / 
               (2*spheres[i].volume));
    spheres[i].semiaxis_c = sqrt(5*(eigenvalue[1] + eigenvalue[0] - eigenvalue[2]) / 
               (2*spheres[i].volume));
    spheres[i].deformation = (spheres[i].semiaxis_a - spheres[i].semiaxis_c) / 
                             (spheres[i].semiaxis_a + spheres[i].semiaxis_c); 
    //double mag = spheres[i].i_eigenvector[0]*spheres[i].i_eigenvector[0] + 
    //             spheres[i].i_eigenvector[1]*spheres[i].i_eigenvector[1] + 
    //             spheres[i].i_eigenvector[2]*spheres[i].i_eigenvector[2];
    //spheres[i].i_eigenvector[0] /= mag;
    //spheres[i].i_eigenvector[1] /= mag;
    //spheres[i].i_eigenvector[2] /= mag;
    double longAxis[3];
    longAxis[0] = spheres[i].i_eigenvector[0];
    longAxis[1] = spheres[i].i_eigenvector[1];
    longAxis[2] = spheres[i].i_eigenvector[2];
    double flowDirec[3];
    flowDirec[0]=1.; flowDirec[1]=0.; flowDirec[2]=0.;
    double cross[3]={0.};
    product(longAxis,flowDirec,cross);
    double crossMag = cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2];
    crossMag = sqrt(crossMag);
    double dot = longAxis[0]*flowDirec[0];
    spheres[i].inclination = atan2(crossMag,dot) * 180.0 / Pi ;
    //printf("D = %f  angle = %f\n",spheres[i].deformation, spheres[i].inclination*180/Pi);
	}
// modified 20160913
//if(spheres[0].volume > 90.51)
//  sphere_pm->V0[0] = 90.51;
//else 
//  sphere_pm->V0[0] = spheres[0].volume + .025;
//if(spheres[0].area > 137.57)
//  sphere_pm->A0[0] = 137.57;
//else 
//  sphere_pm->A0[0] = spheres[0].area + .025;
//
////sphere_pm->V0[0]= spheres[0].volume;
////sphere_pm->A0[0]= spheres[0].area;
//if(spheres[0].avg_sprng_len > 1.)
//  springLength_eq = 1.0;
//else
//  springLength_eq = spheres[0].avg_sprng_len + 0.01;
//
//printf("V0=%f  A0=%f  springLength=%f\n", sphere_pm->V0[0], sphere_pm->A0[0], springLength_eq);

}

void Write_Sphere(struct sphere *spheres, struct monomer *monomers, struct face *faces, 
     struct sphere_param *sphere_pm, int n_step, int addsteps, char *work_dir)
{
	int i, j, k, n, d;
	char filename[160];
	int num_spheres = sphere_pm->Nsphere;
	int num_step;
//	double avg_disp2;
//	double  avg_Rg2;
//	double avg_stretchx, avg_stretchy, avg_stretchz;
//	double avg_dr2;
	double avg_sp_len;
//	double avg_asphr;
//	double avg_acldr;
	static int   buff_nstep[WRITE_BUFFER];
//	static double buff_avgdisp2[WRITE_BUFFER];
//	static double buff_avgdr2[WRITE_BUFFER];
//	static double buff_Rg2[WRITE_BUFFER];
	static double buff_stretch[WRITE_BUFFER][DIMS];
	static double buff_sp_len[WRITE_BUFFER];
//	static double buff_asphr[WRITE_BUFFER];
//	static double buff_acldr[WRITE_BUFFER];
	static double buff_sphereprops[WRITE_BUFFER][10][MAX_DP];
	static double buff_sphereshape[WRITE_BUFFER][11][MAX_DP];
//	static double buff_spherevalues[WRITE_BUFFER][2][MAX_DP];
//	static double buff_spherevolume[WRITE_BUFFER][4][MAX_DP];
	FILE *stream;

	if(n_step == 0) 
  {
		if(addsteps == 0) // what's 'addsteps'? 
		{
			sprintf(filename, "%s/data/avg_props.dat", work_dir);
			stream = fopen(filename, "w");
//			fprintf(stream, "#nstep    disp2    dr2    Rg2    stretch(x,y,z)    sprng_len    asphericity    acylindricity\n");
      fprintf(stream, "#nstep    sprng_len\n");
		  fclose(stream);
			for(n=0; n < num_spheres; n++) 
			{
				sprintf(filename, "%s/data/sphere_props.%d.dat", work_dir, n);
				stream=fopen(filename, "w");
//				fprintf(stream, "#nstep    com    Rg2    avg_sprng_len    disp2    dr2    stretch(x,y,z)    volume    area    theta    omega[z]\n");
				fprintf(stream, "nstep    COM(x,y,z)    volume    area    avgSpringLen    stretch(x,y,z)    omega[z]\n");
        fclose(stream);
				sprintf(filename, "%s/data/sphere_shape.%d.dat", work_dir, n);
				stream=fopen(filename, "w");
//				fprintf(stream, "#nstep    Rx2    Ry2    Rz2    Ia    Ib    Ic\n");
        fprintf(stream, "#nstep    g_eigenX    g_eigenY    g_eigenZ    semiaxis_a    semiaxis_b    semiaxis_c    i_vec_a    deforma   angle\n");
				fclose(stream);
			}
		}
		for(i=0; i < WRITE_BUFFER; i++) 
		{
			//buff_avgdisp2[i] = 0.0;
			//buff_avgdr2[i]   = 0.0;
			//buff_Rg2[i]      = 0.0;
		  buff_sp_len[i]   = 0.0;
			//buf_asphr[i]    = 0.0;
		  //buff_acldr[i]    = 0.0;
			for(d=0; d < DIMS; d++)
				buff_stretch[i][d]=0.0;

			for(k=0; k < MAX_DP; k++)  // Note: MAX_DP==500 
			{
				for(j=0; j < 10; j++)
					buff_sphereprops[i][j][k]=0.0;
				for(j=0; j < 11; j++)
					buff_sphereshape[i][j][k]=0.0;
//				for(j=0; j<2; j++)
//					buff_spherevalues[i][j][k]=0.0;
//				for(j=0; j<4; j++)
//					buff_spherevolume[i][j][k]=0.0;
			}  
		}
	}
	num_step = n_step / sphere_pm->write_time; // buffer index;the variable name is misleading
	if(n_step % sphere_pm->write_time == 0)//Calculate properties and store them in the buffer
	{
	  sphere_props(sphere_pm, spheres, monomers, faces);

//		avg_Rg2      = 0.0;
//		avg_stretchx = 0.0;
//		avg_stretchy = 0.0;
//		avg_stretchz = 0.0;
//		avg_disp2    = 0.0;
//		avg_dr2      = 0.0;
		avg_sp_len   = 0.0;
//		avg_asphr    = 0.0;
//		avg_acldr    = 0.0;

		for(n=0; n < num_spheres; n++) {
//			avg_disp2 += spheres[n].disp2;
//			avg_Rg2   += spheres[n].Rg2;
//			avg_stretchx += spheres[n].stretch[0];
//			avg_stretchy += spheres[n].stretch[1];
//			avg_stretchz += spheres[n].stretch[2];
//			avg_dr2 += spheres[n].dr2;
			avg_sp_len += spheres[n].avg_sprng_len;
//			avg_asphr += spheres[n].asphericity;
//			avg_acldr += spheres[n].acylindricity;
			buff_sphereprops[num_step%WRITE_BUFFER][0][n] = spheres[n].com[0];
			buff_sphereprops[num_step%WRITE_BUFFER][1][n] = spheres[n].com[1];
			buff_sphereprops[num_step%WRITE_BUFFER][2][n] = spheres[n].com[2];
      buff_sphereprops[num_step%WRITE_BUFFER][3][n] = spheres[n].volume;
			buff_sphereprops[num_step%WRITE_BUFFER][4][n] = spheres[n].area;
      buff_sphereprops[num_step%WRITE_BUFFER][5][n] = spheres[n].avg_sprng_len;
      buff_sphereprops[num_step%WRITE_BUFFER][6][n] = spheres[n].stretch[0];
			buff_sphereprops[num_step%WRITE_BUFFER][7][n] = spheres[n].stretch[1];
			buff_sphereprops[num_step%WRITE_BUFFER][8][n] = spheres[n].stretch[2];
      buff_sphereprops[num_step%WRITE_BUFFER][9][n] = spheres[n].omega[2];
//			buff_sphereprops[num_step%WRITE_BUFFER][3][n] = spheres[n].Rg2;
//			buff_sphereprops[num_step%WRITE_BUFFER][4][n] = spheres[n].avg_sprng_len;
//			buff_sphereprops[num_step%WRITE_BUFFER][5][n] = spheres[n].disp2;
//			buff_sphereprops[num_step%WRITE_BUFFER][6][n] = spheres[n].dr2;
//			buff_sphereprops[num_step%WRITE_BUFFER][7][n] = spheres[n].stretch[0];
//			buff_sphereprops[num_step%WRITE_BUFFER][8][n] = spheres[n].stretch[1];
//			buff_sphereprops[num_step%WRITE_BUFFER][9][n] = spheres[n].stretch[2];
//			buff_sphereprops[num_step%WRITE_BUFFER][10][n] = spheres[n].volume;
//			buff_sphereprops[num_step%WRITE_BUFFER][11][n] = spheres[n].area;
//			buff_sphereprops[num_step%WRITE_BUFFER][12][n] = spheres[n].theta;
//			buff_sphereprops[num_step%WRITE_BUFFER][13][n] = spheres[n].omega[2];
//			buff_sphereshape[num_step%WRITE_BUFFER][0][n] = spheres[n].Rx2;
//			buff_sphereshape[num_step%WRITE_BUFFER][1][n] = spheres[n].Ry2;
//			buff_sphereshape[num_step%WRITE_BUFFER][2][n] = spheres[n].Rz2;
//			buff_sphereshape[num_step%WRITE_BUFFER][3][n] = spheres[n].Ia;
//			buff_sphereshape[num_step%WRITE_BUFFER][4][n] = spheres[n].Ib;
//			buff_sphereshape[num_step%WRITE_BUFFER][5][n] = spheres[n].Ic;
			buff_sphereshape[num_step%WRITE_BUFFER][0][n] = spheres[n].g_eigenvalue[0];
			buff_sphereshape[num_step%WRITE_BUFFER][1][n] = spheres[n].g_eigenvalue[1];
			buff_sphereshape[num_step%WRITE_BUFFER][2][n] = spheres[n].g_eigenvalue[2];
			buff_sphereshape[num_step%WRITE_BUFFER][3][n] = spheres[n].semiaxis_a;
			buff_sphereshape[num_step%WRITE_BUFFER][4][n] = spheres[n].semiaxis_b;
			buff_sphereshape[num_step%WRITE_BUFFER][5][n] = spheres[n].semiaxis_c;
			buff_sphereshape[num_step%WRITE_BUFFER][6][n] = spheres[n].i_eigenvector[0];
			buff_sphereshape[num_step%WRITE_BUFFER][7][n] = spheres[n].i_eigenvector[1];
			buff_sphereshape[num_step%WRITE_BUFFER][8][n] = spheres[n].i_eigenvector[2];
			buff_sphereshape[num_step%WRITE_BUFFER][9][n] = spheres[n].deformation;
			buff_sphereshape[num_step%WRITE_BUFFER][10][n] = spheres[n].inclination;
			//buff_sphereshape[num_step%WRITE_BUFFER][11][n] = spheres[n].i_eigenvector[5];
			//buff_sphereshape[num_step%WRITE_BUFFER][12][n] = spheres[n].i_eigenvector[6];
			//buff_sphereshape[num_step%WRITE_BUFFER][13][n] = spheres[n].i_eigenvector[7];
			//buff_sphereshape[num_step%WRITE_BUFFER][14][n] = spheres[n].i_eigenvector[8];
		}
//		avg_disp2    /= num_spheres;
//		avg_Rg2      /= num_spheres;
//		avg_stretchx /= num_spheres;
//		avg_stretchy /= num_spheres;
//		avg_stretchz /= num_spheres;
//		avg_dr2      /= num_spheres;
		avg_sp_len   /= num_spheres;
//		avg_asphr    /= num_spheres;
//		avg_acldr    /= num_spheres;

		buff_nstep[num_step % WRITE_BUFFER] = (int)((n_step+addsteps)*sphere_pm->MD_steps*
                                          sphere_pm->dt);
//		buff_avgdisp2[num_step % WRITE_BUFFER] = avg_disp2;
//		buff_avgdr2[num_step%WRITE_BUFFER]     = avg_dr2;
//		buff_Rg2[num_step%WRITE_BUFFER]        = avg_Rg2;
//		buff_stretch[num_step%WRITE_BUFFER][0] = avg_stretchx;
//		buff_stretch[num_step%WRITE_BUFFER][1] = avg_stretchy;
//		buff_stretch[num_step%WRITE_BUFFER][2] = avg_stretchz;
		buff_sp_len[num_step%WRITE_BUFFER]     = avg_sp_len;
//		buff_asphr[num_step%WRITE_BUFFER]      = avg_asphr;
//		buff_acldr[num_step%WRITE_BUFFER]      = avg_acldr;

		if((num_step+1) % WRITE_BUFFER == 0)  // Write out all the data storing in the buffer.
		{
			sprintf(filename, "%s/data/avg_props.dat", work_dir);
			stream = fopen(filename, "a");
			for(i=0; i < WRITE_BUFFER; i++)
			{
//				fprintf(stream, "%d %le %le %le %le %le %le %le %le %le\n", 
//						buff_nstep[i], buff_avgdisp2[i], buff_avgdr2[i], buff_Rg2[i], 
//						buff_stretch[i][0], buff_stretch[i][1], buff_stretch[i][2], 
//						buff_sp_len[i], buff_asphr[i], buff_acldr[i]);
				fprintf(stream, "%d %le\n",buff_nstep[i], buff_sp_len[i]);	
			}
			fclose(stream);
			for(n=0; n < num_spheres; n++) {	
				for(j=0; j < WRITE_BUFFER; j++) 
        {
					sprintf(filename, "%s/data/sphere_props.%d.dat", work_dir, n);
					stream=fopen(filename, "a");
					fprintf(stream, "%d ", buff_nstep[j]);
//					for(i=0; i<14; i++)
//						fprintf(stream, "%le ", buff_sphereprops[j][i][n]);
          for(i=0; i < 10; i++)
            fprintf(stream, "%le ", buff_sphereprops[j][i][n]);
					fprintf(stream, "\n");
					fclose(stream);
					sprintf(filename, "%s/data/sphere_shape.%d.dat", work_dir, n);
					stream=fopen(filename, "a");
					fprintf(stream, "%d ", buff_nstep[j]);
					for(i=0; i < 11; i++)
						fprintf(stream, "%le ", buff_sphereshape[j][i][n]);
					fprintf(stream, "\n");
					fclose(stream);
				}
			}
		}
	}
}

void Write_Fluid(Float ***velcs_df, int **node_map, int nstep, char *work_dir)
{
	extern int max_x, max_y, max_z;
	extern int num_x;
	char filename[160];
	double prop[DIMS+1+2*DIMS];  // properties for fluid density, 
	double lvx, lvy, lvz;        // fluid velocity, and fluid stress.
	FILE *stream;
	long i,j,k, p,q, q1, nxy;

	/* Write out velocity field in vtk format */
	write_velocity_field(nstep, velcs_df, node_map, work_dir);  
}

/* Output velocity field in the xy-plane in vtk format */
void write_velocity_field(int istep,Float ***velcs_df,int **node_map,char *work_dir)
{
	char PREFIX[20]="u";
	char FILE_NAME[100];
	FILE *outfile;
	int i,j,k,l,d;
	int q, nxy;
	double v[3];
	extern int max_x;
	extern int max_y;
	extern int max_z;

	if(istep < 10)
		sprintf(FILE_NAME, "%s/data/%s.000%d.vtk", work_dir, PREFIX, istep);
	else if(istep < 100)
		sprintf(FILE_NAME, "%s/data/%s.00%d.vtk", work_dir, PREFIX, istep);
	else if(istep < 1000)
		sprintf(FILE_NAME, "%s/data/%s.0%d.vtk", work_dir, PREFIX, istep);
	else
		sprintf(FILE_NAME, "%s/data/%s.%d.vtk", work_dir, PREFIX, istep);

	//sprintf(FILE_NAME, "%s/data/u%d.dat", work_dir, istep);
	outfile = fopen(FILE_NAME, "w");
	fprintf(outfile, "# vtk DataFile Version 2.3   \n");
	fprintf(outfile, "title u field(istep=%d)      \n", istep);
	fprintf(outfile, "ASCII                        \n");
	fprintf(outfile, "DATASET STRUCTURED_POINTS    \n");
	//fprintf(outfile, "DIMENSIONS %d %d %d\n",max_x, max_y, 1);
  fprintf(outfile, "DIMENSIONS %d %d %d \n", max_x, max_y, max_z);
	fprintf(outfile, "ORIGIN  %le %le %le\n", 0.0, 0.0, 0.0);
	fprintf(outfile, "SPACING %le %le %le\n", 1.0, 1.0, 1.0);
	//fprintf(outfile, "POINT_DATA %d\n", (max_x)*max_y);
  fprintf(outfile, "POINT_DATA %d\n", max_x*max_y*max_z);
	fprintf(outfile, "VECTORS outerFluid double      \n"); // 20161228
//  fprintf(outfile, "VECTORS fluid double      \n"); //20170120
//	fclose(outfile);
//	outfile = fopen(FILE_NAME, "a");
	for(k=1; k <= max_z; k++) 
  {
    //fprintf(outfile, "VECTORS velSlice(z=%d) double      \n", k); // 20161228
		for(j=1; j <= max_y; j++) {
			for(i=1; i <= max_x; i++) {	
				nxy = i * (max_y+2) + j;
				v[0] = 0.0;
				v[1] = 0.0;
				v[2] = 0.0;
				if(node_map[nxy][k]==0) { // Fluid node
					for(q=0; q < Num_Dir; q++) {
						v[0] += c_x[q]*velcs_df[nxy][q][k];
						v[1] += c_y[q]*velcs_df[nxy][q][k];
						v[2] += c_z[q]*velcs_df[nxy][q][k];
					}
				}  //20170120
        //for(d=0; d<3; d++)
        //  v[d] /= 36.0;  // divided by fluid density; Fluid density is defined in drive.c.

				//	fprintf(outfile, "(%- d\t %- d\t %- d)\t%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", i,j,k,v[0], v[1], v[2]); 
				fprintf(outfile, "%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", v[0], v[1], v[2]);
				//	fprintf(outfile, "%- d\t%- d\t%- d\n", i, j, k); 
			}
		}
	}
//20170120

	fprintf(outfile, "VECTORS innerFluid double      \n"); // 20161227
	for(k=1; k <= max_z; k++) 
  {
    //fprintf(outfile, "VECTORS innerFluid(z=%d) double      \n",k);
		for(j=1; j <= max_y; j++) {
			for(i=1; i <= max_x; i++) {	
				nxy = i * (max_y+2) + j;
				v[0] = 0.0;
				v[1] = 0.0;
				v[2] = 0.0;
				if(node_map[nxy][k]==2) { // inner fluid node
					for(q=0; q < Num_Dir; q++) {
						v[0] += c_x[q]*velcs_df[nxy][q][k];
						v[1] += c_y[q]*velcs_df[nxy][q][k];
						v[2] += c_z[q]*velcs_df[nxy][q][k];
					}
				}
        //for(d=0; d<3; d++)
        //  v[d] /= 36.0;  // divided by fluid density; fluid density is defined in drive.c
				fprintf(outfile, "%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", v[0], v[1], v[2]);
      }
		}
	}

	fclose(outfile);

	/*   if(istep < 10) */
	/*     sprintf(FILE_NAME, "%s/%s.000%d.dat", work_dir, PREFIX, istep); */
	/*   else if(istep < 100) */
	/*     sprintf(FILE_NAME, "%s/%s.00%d.dat", work_dir, PREFIX, istep); */
	/*   else if(istep < 1000) */
	/*     sprintf(FILE_NAME, "%s/%s.0%d.dat", work_dir, PREFIX, istep); */
	/*   else */
	/*     sprintf(FILE_NAME, "%s/%s.%d.dat", work_dir, PREFIX, istep); */

	/*   outfile = fopen(FILE_NAME, "w"); */

	/*   fprintf(outfile," VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"W\"\n"); */
	/*   fprintf(outfile, "ZONE I=%ld, J=%ld\n", max_x, max_y); */
	/*   k=max_z/2; */
	/*   for(j=1; j <= max_y; j++) { */
	/*     for(i=1; i<= max_x; i++) { */
	/* 	nxy = i*(max_y+2) + j; */

	/* 	if (node_map[nxy][k]==0) {         /\* Fluid node *\/ */
	/* 	  v[0]=0.0; */
	/* 	  v[1]=0.0; */
	/* 	  v[2]=0.0; */
	/* 	  for(q=0; q < Num_Dir; q++) { */
	/* 	    v[0]+= c_x[q]*velcs_df[nxy][q][k]; */
	/* 	    v[1]+= c_y[q]*velcs_df[nxy][q][k]; */
	/* 	    v[2]+= c_z[q]*velcs_df[nxy][q][k]; */
	/* 	  } */

	/* 	  for(d=0; d<DIMS; d++) */
	/* 	    v[d]/=Rho_Fl; */

	/* 	  fprintf(outfile, "%d %d %le %le %le \n", i,j, v[0], v[1], v[2]); */
	/* 	} */
	/*     } */
	/*   } */
	/*   fclose(outfile); */

	/* write out system fluid density and momentum */

	/*   if(nstep % WRITE_BUFFER == 0){ */
	/*     sprintf(filename, "%s/data/stress.dat", work_dir); */
	/*     stream = fopen(filename, "a"); */

	/*     i = (max_x+1)/2; */
	/*     j = (max_y+1)/2; */
	/*     k = (max_z+1)/2; */

	/*     for(p=0; p<DIMS+1+2*DIMS; p++) */
	/*       prop[p]=0.0; */

	/*     for(i=1; i<=num_x; i++) { */
	/*       for(j=0; j<max_y; j++) { */
	/* 	nxy = i*max_y + j; */

	/* 	for(k=0; k<max_z; k++) { */

	/* 	  lvx =0.0; */
	/* 	  lvy =0.0; */
	/* 	  lvz =0.0; */
	/* 	  for(q=0; q < Num_Dir; q++) { */
	/* 	  /\*	    prop[0] += evector[0][q]*velcs_df[nxy][q][k]; *\/ */
	/* 	    lvx += evector[1][q]*velcs_df[nxy][q][k]; */
	/* 	    lvy += evector[2][q]*velcs_df[nxy][q][k]; */
	/* 	    lvz += evector[3][q]*velcs_df[nxy][q][k]; */

	/* 	    /\* stresses *\/ */
	/* 	    prop[4] += evector[4][q]*velcs_df[nxy][q][k]; */
	/* 	    prop[5] += evector[5][q]*velcs_df[nxy][q][k]; */
	/* 	    prop[6] += evector[6][q]*velcs_df[nxy][q][k]; */
	/* 	    prop[7] += evector[7][q]*velcs_df[nxy][q][k]; */
	/* 	    prop[8] += evector[8][q]*velcs_df[nxy][q][k]; */
	/* 	    prop[9] += evector[9][q]*velcs_df[nxy][q][k]; */

	/* 	  } */

	/* 	  prop[0] += lvx*lvx+lvy*lvy+lvz*lvz; */
	/* 	  prop[1] += lvx; */
	/* 	  prop[2] += lvy; */
	/* 	  prop[3] += lvz; */

	/* 	} */
	/*       } */
	/*     } */


	/*     for(p=0; p<=DIMS; p++) */
	/*       prop[p] /= (num_x*max_y*max_z); */


	/*     fprintf(stream, "%d %le %le %le %le ", nstep, prop[0], prop[1], prop[2], prop[3]); */
	/*     fprintf(stream, "%le %le %le %le %le %le\n", prop[4], prop[5], prop[6], prop[7], prop[8], prop[9]); */

	/*     fclose(stream); */
	/*   } */

}

void Write_Velslice(struct monomer *monomers, Float ***velcs_df, int **node_map, 
     int nstep, char *work_dir)
{
	/* Write out the velocity slicing through the monomer position */
	extern int max_x, max_y, max_z;
	int i,j,k,p,q, nxy;
	double lvx, lvy, lvz;
	double prop[DIMS+1];
	char filename[160];
	FILE *stream;
	if(nstep < 10)
		sprintf(filename, "%s/data/velslice.000%d.dat", work_dir, nstep);
	else if(nstep < 100)
		sprintf(filename, "%s/data/velslice.00%d.dat", work_dir, nstep);
	else if(nstep < 1000)
		sprintf(filename, "%s/data/velslice.0%d.dat", work_dir, nstep);
	else if(nstep >= 1000)
		sprintf(filename, "%s/data/velslice.%d.dat", work_dir, nstep);

	stream = fopen(filename, "w");

	//i = floor(monomers[0].pos[0]); 
	//k = floor(monomers[0].pos[2]); 

	i = (max_x+1)/2;
	k = (max_z+1)/2;
	// j = (max_y+1)/2;

	fprintf(stream, "#monpos (%d, j, %d)\n", i, k);
	for(j=1; j <= max_y; j++) 
	{  
		nxy = i*(max_y+2) + j;

		for(p=0; p<DIMS+1; p++)
			prop[p]=0.0;

		if(node_map[nxy][k]==0) // Fluid node
		{         
			lvx =0.0;
			lvy =0.0;
			lvz =0.0;
			for(q=0; q < Num_Dir; q++) 
			{
				prop[0] += evector[0][q]*velcs_df[nxy][q][k];
				lvx += c_x[q]*velcs_df[nxy][q][k];
				lvy += c_y[q]*velcs_df[nxy][q][k];
				lvz += c_z[q]*velcs_df[nxy][q][k];
			}

			/*prop[0] += lvx*lvx+lvy*lvy+lvz*lvz; */
			prop[1] = lvx /Rho_Fl;
			prop[2] = lvy /Rho_Fl;
			prop[3] = lvz /Rho_Fl;
			fprintf(stream, "%d %le %le %le %le \n", j, prop[0], prop[1], prop[2], prop[3]);
		}
	}
	fclose(stream);     
}

/* Output velocity field in the xy-plane in vtk format */
void Write_Nodemap(int **node_map, int istep, char *work_dir)
{
	char PREFIX[20]="node_map";
	char FILE_NAME[100];
	FILE *outfile;
	int i,j,k,l,d;
	int q, nxy;
	double v[3];
	extern int max_x;
	extern int max_y;
	extern int max_z;

	if(istep < 10)
		sprintf(FILE_NAME, "%s/data/%s.000%d.vtk", work_dir, PREFIX, istep);
	else if(istep < 100)
		sprintf(FILE_NAME, "%s/data/%s.00%d.vtk", work_dir, PREFIX, istep);
	else if(istep < 1000)
		sprintf(FILE_NAME, "%s/data/%s.0%d.vtk", work_dir, PREFIX, istep);
	else
		sprintf(FILE_NAME, "%s/data/%s.%d.vtk", work_dir, PREFIX, istep);

	outfile = fopen(FILE_NAME, "w");
	fprintf(outfile, "# vtk DataFile Version 2.3   \n");
	fprintf(outfile, "title u field(istep=%d)      \n", istep);
	fprintf(outfile, "ASCII                        \n");
	fprintf(outfile, "DATASET STRUCTURED_POINTS    \n");
	fprintf(outfile, "DIMENSIONS %d %d %d \n", max_x, max_y, max_z);
	fprintf(outfile, "ORIGIN  %le %le %le\n", 0.0, 0.0, 0.0);
	fprintf(outfile, "SPACING %le %le %le\n", 1.0, 1.0, 1.0);
	fprintf(outfile, "POINT_DATA %d\n", max_x*max_y*max_z);
	fprintf(outfile, "VECTORS velocity double      \n");
	fclose(outfile);
	outfile = fopen(FILE_NAME, "a");
	for(k=1; k<=max_z; k++) 
	{
		for(j=1; j<=max_y; j++)
			for(i=1; i<=max_x; i++) 
			{	
				nxy = i*(max_y+2)+j;
				v[0] = 0.0;
				v[1] = 0.0;
				v[2] = 0.0;
				if (node_map[nxy][k]==2)  /* Fluid node */ 
				{          
					v[0] = 0.1;
					v[1] = 0.1;
					v[2] = 0.1;
				}
				/*	fprintf(outfile, "(%- d\t %- d\t %- d)\t%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", i,j,k,v[0], v[1], v[2]); */
				fprintf(outfile, "%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", v[0], v[1], v[2]);
				/*	fprintf(outfile, "%- d\t%- d\t%- d\n", i, j, k); */
			}
	}
	fclose(outfile);
}


//void Write_Monomer(struct sphere *spheres, struct monomer *monomers, struct sphere_param 
//		*sphere_pm, int n_step, int addsteps, char *work_dir)
//{
//	extern int max_x, max_y, max_z;
//	static int buff_nstep[WRITE_BUFFER];
//	static double buff_avgdisp2[WRITE_BUFFER], buff_avgdr2[WRITE_BUFFER], buff_avgvel[WRITE_BUFFER];
//	static double buff_monvel[WRITE_BUFFER][500], buff_tempscale[WRITE_BUFFER];
//	char filename[160];
//	int j;
//	int n, n2, d;
//	int maxsize[DIMS];
//	int num_write, num_step;
//	int tot_time, add_time;
//	double avg_disp2, avg_vel, avg_E, avg_dr2;
//	int num_beads = sphere_pm->num_beads;
//	FILE *stream;
//
//	maxsize[0]=max_x;
//	maxsize[1]=max_y;
//	maxsize[2]=max_z;
//	tot_time = (n_step+addsteps)*sphere_pm->MD_steps*sphere_pm->dt;
//	add_time = addsteps*sphere_pm->MD_steps*sphere_pm->dt;
//
//	if(n_step == 0) 
//	{
//		for(n=0; n<num_beads; n++)
//		{
//			for(d=0; d<DIMS; d++) 
//			{
//				monomers[n].pos0[d]    = monomers[n].pos_pbc[d];
//				monomers[n].pos_old[d] = monomers[n].pos_pbc[d];
//			}
//		}
//		for(j=0; j<WRITE_BUFFER; j++) 
//		{
//			buff_nstep[j]=0;
//			buff_avgdisp2[j]=0.0;
//			buff_avgdr2[j]=0.0;
//			buff_avgvel[j]=0.0;
//			buff_tempscale[j]=0.0;
//		}
//	}
//
//	if(sphere_pm->numsteps > 1001) 
//		num_write = sphere_pm->relax_time/(sphere_pm->MD_steps*sphere_pm->dt);
//	else if(sphere_pm->num_beads == 1) 
//		num_write = 10;
//	else 
//		num_write = 1;
//
//	if((n_step+addsteps) % num_write == 0) 
//	{
//		/* Write out monomer position and velocity */
//		sprintf(filename, "%s/data/monpos.dat", work_dir);
//		stream=fopen(filename, "a");
//
//		if(num_beads < 500)
//		{ 
//			for(n=0; n<num_beads; n++)
//			{  
//				fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\t %12.10le %12.10le %12.10le\n", 
//						tot_time, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2], 
//						monomers[n].vel[0], monomers[n].vel[1], monomers[n].vel[2], 
//						monomers[n].fluid_vel[0], monomers[n].fluid_vel[1], monomers[n].fluid_vel[2]);
//			}
//		}
//		else
//		{
//			for(n=0; n<num_beads; n++)
//			{ 
//				fprintf(stream, "%d %12.10le %12.10le %12.10le\n", 
//						tot_time, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2]);
//			}
//		}  
//		fclose(stream);
//
//		/* write out bonding config */
//		if (WRITE_BOND)               // By default (head.h), WRITE_BOND==0!
//		{
//			sprintf(filename, "%s/data/bond.dat", work_dir);
//			stream=fopen(filename, "a");
//
//			for (n=0 ; n < num_beads ; n++)
//			{  
//				for(j=1 ; j <= monomers[n].blist[0][0] ; j++) 
//				{
//					n2 = monomers[n].blist[j][0];
//
//					fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\n", 
//							tot_time, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2], 
//							monomers[n2].pos[0], monomers[n2].pos[1], monomers[n2].pos[2]);
//				}
//			}
//
//			fclose(stream);
//		}
//	}
//
//	/* write out monomer velocities */
//	if(num_beads < 5)
//		for(n=0; n<num_beads; n++) {
//			sprintf(filename, "%s/data/monvel_%d.dat", work_dir, n);
//			stream=fopen(filename, "a");
//
//			fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le\n", tot_time, monomers[n].vel[0], monomers[n].vel[1], monomers[n].vel[2], monomers[n].fricforce[0]);
//
//			fclose(stream);
//		}
//
//	num_step = n_step / sphere_pm->write_time;
//	if(n_step % sphere_pm->write_time == 0) {
//		avg_disp2 = 0.0;
//		avg_vel = 0.0;
//		avg_dr2 = 0.0;
//
//		for(n=0; n < num_beads; n++) {
//			for(d=0; d < DIMS; d++) {
//				avg_disp2 += (monomers[n].pos_pbc[d] - monomers[n].pos0[d]) * 
//                     (monomers[n].pos_pbc[d] - monomers[n].pos0[d]);
//				avg_vel += monomers[n].vel[d] * monomers[n].vel[d];
//
//				/* Calculate monomer mean square displacement */
//				monomers[n].dr2 += (monomers[n].pos_pbc[d]-monomers[n].pos_old[d])*
//                           (monomers[n].pos_pbc[d]-monomers[n].pos_old[d]);
//				monomers[n].pos_old[d] = monomers[n].pos_pbc[d];
//			}
//			avg_dr2 += monomers[n].dr2;
//		}
//		avg_disp2 /= num_beads;
//		avg_vel /= num_beads;
//		avg_dr2 /= num_beads;
//
//		buff_nstep[num_step%WRITE_BUFFER]=tot_time;
//		buff_avgdisp2[num_step%WRITE_BUFFER]=avg_disp2;
//		buff_avgdr2[num_step%WRITE_BUFFER]=avg_dr2;
//		buff_avgvel[num_step%WRITE_BUFFER]=avg_vel*Rho_Fl*sphere_pm->monmass;  
//		buff_tempscale[num_step%WRITE_BUFFER]=sphere_pm->tempscale;  
//
//		if((num_step+1)%WRITE_BUFFER == 0) {
//			sprintf(filename, "%s/data/avg_disp2.dat", work_dir);
//			stream = fopen(filename, "a");
//			for(j=0; j<WRITE_BUFFER; j++)
//				fprintf(stream, "%d %le %le %le\n", buff_nstep[j], buff_avgdisp2[j], buff_avgdr2[j], buff_avgvel[j]);
//			fclose(stream);
//		}
//	} 
//}

void Write_Sphereconfig(struct sphere *spheres, struct monomer *monomers, struct face 
		*faces, struct sphere_param *sphere_pm, int n_step, int addsteps,	char *work_dir)
{
	extern int max_x, max_y, max_z;
	int i,j,k,type;
	int nbead[NTYPES], nface[NTYPES];
	int bead_0, bead_f, face_0, face_f;
	int nface_temp;
	int spherenum;
	double maxsize[DIMS];
	char filename[160];
	int n = (addsteps+n_step)/sphere_pm->write_config;
	int m = (addsteps+n_step)%sphere_pm->write_config;
	FILE *outfile;
	maxsize[0] = max_x;
	maxsize[1] = max_y;
	maxsize[2] = max_z;

	if(m == 0) 
	{
		for(type = 0; type < 1 /*NTYPES*/; type++) // Modification 20170714: Temporary  
		{                                          
			nbead[type] = sphere_pm->N_per_sphere[type]*sphere_pm->Ntype[type];
			nface[type] = sphere_pm->face_per_sphere[type]*sphere_pm->Ntype[type];

			if(type ==0) 
			{
				bead_0 = 0;
				face_0 = 0;

				bead_f = nbead[0];
				face_f = nface[0];
			}
			else 
			{
				bead_0 = nbead[type-1];
				face_0 = nface[type-1];

				bead_f += nbead[type];
				face_f += nface[type];
			}

			if(n < 10)
				sprintf(filename, "%s/data/bond%d_t000%d.vtk", work_dir, type, n);
			else if(n < 100)
				sprintf(filename, "%s/data/bond%d_t00%d.vtk", work_dir, type, n);
			else if(n < 1000)
				sprintf(filename, "%s/data/bond%d_t0%d.vtk", work_dir, type, n);
			else if(n >= 1000)
				sprintf(filename, "%s/data/bond%d_t%d.vtk", work_dir, type, n);

			outfile = fopen(filename, "w");
			fprintf(outfile, "# vtk DataFile Version 2.3   \n");
			fprintf(outfile, "title red blood cell data %d \n", n);
			fprintf(outfile, "ASCII                        \n\n");
			fprintf(outfile, "DATASET UNSTRUCTURED_GRID \n");
			fprintf(outfile, "POINTS %d float\n", nbead[type]);

			for(i=bead_0; i<bead_f; i++)
				fprintf(outfile, "%8.6lf %8.6lf %8.6lf\n", monomers[i].pos[0],monomers[i].pos[1],monomers[i].pos[2]);

			/* check if the face spans across periodic boundaries */
			nface_temp = 0;
			for(i=face_0; i<face_f; i++) {
				if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
					if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
						if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
							if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
								if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
									if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
										if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
											if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
												if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
													nface_temp++;
			}
			fprintf(outfile, "CELLS %d %d\n", nface_temp, nface_temp*4);

			for(i=face_0; i<face_f; i++) 
				if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
					if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
						if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
							if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
								if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
									if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
										if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
											if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
												if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
													fprintf(outfile, "3 %d %d %d\n", faces[i].vertices[0]-bead_0, faces[i].vertices[1]-bead_0, 
                                  faces[i].vertices[2]-bead_0);

			fprintf(outfile, "CELL_TYPES %d\n", nface_temp);
			for(i=0; i<nface_temp; i++)
				fprintf(outfile, "5\n");
			fclose(outfile);
		}
	}
}

void write_growth_config(struct sphere *spheres, struct monomer *monomers, struct face 
		 *faces, struct sphere_param *sphere_pm, int n_step, int addsteps,	char *work_dir)
{
	extern int max_x, max_y, max_z;
	int i,j,k,type;
	int nbead[NTYPES], nface[NTYPES];
	int bead_0, bead_f, face_0, face_f;
	int nface_temp;
	int spherenum;
	double maxsize[DIMS];
	char filename[160];
	int n = (addsteps+n_step)/sphere_pm->write_config;
	int m = (addsteps+n_step)%sphere_pm->write_config;
	FILE *outfile;
	maxsize[0] = max_x;
	maxsize[1] = max_y;
	maxsize[2] = max_z;

	if(m == 0) 
	{
		for(type = 0; type < 1 /*NTYPES*/; type++)  // Modification 20170714: Temporary
		{                                          
			nbead[type] = sphere_pm->N_per_sphere[type]*sphere_pm->Ntype[type];
			nface[type] = sphere_pm->face_per_sphere[type]*sphere_pm->Ntype[type];

			if(type ==0) 
			{
				bead_0 = 0;
				face_0 = 0;

				bead_f = nbead[0];
				face_f = nface[0];
			}
			else 
			{
				bead_0 = nbead[type-1];
				face_0 = nface[type-1];

				bead_f += nbead[type];
				face_f += nface[type];
			}

			if(n < 10)
				sprintf(filename, "%s/data/growth%d_t000%d.vtk", work_dir, type, n);
			else if(n < 100)
				sprintf(filename, "%s/data/growth%d_t00%d.vtk", work_dir, type, n);
			else if(n < 1000)
				sprintf(filename, "%s/data/growth%d_t0%d.vtk", work_dir, type, n);
			else if(n >= 1000)
				sprintf(filename, "%s/data/growth%d_t%d.vtk", work_dir, type, n);

			outfile = fopen(filename, "w");
			fprintf(outfile, "# vtk DataFile Version 2.3   \n");
			fprintf(outfile, "title red blood cell data %d \n", n);
			fprintf(outfile, "ASCII                        \n\n");
			fprintf(outfile, "DATASET UNSTRUCTURED_GRID \n");
			fprintf(outfile, "POINTS %d float\n", nbead[type]);

			for(i=bead_0; i<bead_f; i++)
				fprintf(outfile, "%8.6lf %8.6lf %8.6lf\n", monomers[i].pos[0],monomers[i].pos[1],monomers[i].pos[2]);

			/* check if the face spans across periodic boundaries */
			nface_temp = 0;
			for(i=face_0; i<face_f; i++) {
				if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
					if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
						if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
							if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
								if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
									if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
										if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
											if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
												if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
													nface_temp++;
			}
			fprintf(outfile, "CELLS %d %d\n", nface_temp, nface_temp*4);

			for(i=face_0; i<face_f; i++) 
				if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
					if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
						if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
							if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
								if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
									if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
										if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
											if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
												if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
													fprintf(outfile, "3 %d %d %d\n", faces[i].vertices[0]-bead_0, faces[i].vertices[1]-bead_0, 
                                  faces[i].vertices[2]-bead_0);

			fprintf(outfile, "CELL_TYPES %d\n", nface_temp);
			for(i=0; i<nface_temp; i++)
			  fprintf(outfile, "5\n");
 			fclose(outfile);
		}
	}
}

void Write_overlap_config (int n, struct monomer *monomers, struct face *faces, struct sphere_param *sphere_pm, char *work_dir)
{
	extern int max_x, max_y, max_z;
	int i,j,k,type;
	int nbead[NTYPES], nface[NTYPES];
	int bead_0, bead_f, face_0, face_f;
	int nface_temp;
	int spherenum;
	double maxsize[DIMS];
	char filename[160];
	//int n = (addsteps+n_step)/sphere_pm->write_config;
	//int m = (addsteps+n_step)%sphere_pm->write_config;
	FILE *outfile;
	maxsize[0] = max_x;
	maxsize[1] = max_y;
	maxsize[2] = max_z;

	//if(m == 0) 
	//{
		for(type = 0; type < 1/*NTYPES*/; type++)   // only single type of particles
		{                                          
			nbead[type] = sphere_pm->N_per_sphere[type]*sphere_pm->Ntype[type];
			nface[type] = sphere_pm->face_per_sphere[type]*sphere_pm->Ntype[type];
			if(type ==0) {
				bead_0 = 0;
				face_0 = 0;
				bead_f = nbead[0];
				face_f = nface[0];
			}
			else {
				bead_0 = nbead[type-1];
				face_0 = nface[type-1];
				bead_f += nbead[type];
				face_f += nface[type];
			}

			if(n < 10)
				sprintf(filename, "%s/data/overlap%d_t000%d.vtk", work_dir, type, n);
			else if(n < 100)
				sprintf(filename, "%s/data/overlap%d_t00%d.vtk", work_dir, type, n);
			else if(n < 1000)
				sprintf(filename, "%s/data/overlap%d_t0%d.vtk", work_dir, type, n);
			else if(n >= 1000)
				sprintf(filename, "%s/data/overlap%d_t%d.vtk", work_dir, type, n);

			outfile = fopen(filename, "w");
			fprintf(outfile, "# vtk DataFile Version 2.3   \n");
			fprintf(outfile, "title red blood cell data %d \n", n);
			fprintf(outfile, "ASCII                        \n\n");
			fprintf(outfile, "DATASET UNSTRUCTURED_GRID \n");
			fprintf(outfile, "POINTS %d float\n", nbead[type]);

			for(i = bead_0; i < bead_f; i++)
				fprintf(outfile, "%8.6lf %8.6lf %8.6lf\n", monomers[i].pos_pbc[0],monomers[i].pos_pbc[1],monomers[i].pos_pbc[2]);

			/* check if the face spans across periodic boundaries */
			nface_temp = 0;
			for(i = face_0; i < face_f; i++) {
				if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
					if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
						if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
							if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
								if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
									if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
										if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
											if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
												if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
													nface_temp++;
			}
			fprintf(outfile, "CELLS %d %d\n", nface_temp, nface_temp*4);

			for(i = face_0; i < face_f; i++) 
				if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
					if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
						if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
							if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
								if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
									if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
										if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
											if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
												if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
													fprintf(outfile, "3 %d %d %d\n", faces[i].vertices[0]-bead_0, faces[i].vertices[1]-bead_0, 
                                  faces[i].vertices[2]-bead_0);
			fprintf(outfile, "CELL_TYPES %d\n", nface_temp);
			for(i=0; i < nface_temp; i++)
				fprintf(outfile, "5\n");

      fprintf(outfile,"POINT_DATA %d\n",nbead[type]);
      fprintf(outfile,"VECTORS vertexFluVel double\n");
      for(i = bead_0; i < bead_f; i++)
//        fprintf(outfile, "%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", monomers[i].fluid_vel[0], 
//                monomers[i].fluid_vel[1], monomers[i].fluid_vel[2]);
        fprintf(outfile, "%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", monomers[i].v_diff[0], 
                monomers[i].v_diff[1], monomers[i].v_diff[2]); //20170125
			fclose(outfile);
		}
	//}
}



void check_vel_conservation(Float ***velcs_df, int nstep)
{
	extern int max_x, max_y, max_z;
	int i, j, k, d, q;
	int nxy;
	Float total_deltaf[DIMS];

	i=(1+max_x)/2;
	j=(1+max_y)/2;
	k=(1+max_z)/2;

	for(d=0; d<DIMS; d++) 
		total_deltaf[d]=0.0;

	nxy = i*max_y+j;

	for(d=0; d<DIMS; d++) {
		for(q=0; q<Num_Dir; q++) {
			total_deltaf[0] +=velcs_df[nxy][q][k]*c_x[q];
			total_deltaf[1] +=velcs_df[nxy][q][k]*c_y[q];
			total_deltaf[2] +=velcs_df[nxy][q][k]*c_z[q];
		}
	}

	for(d=0; d<DIMS; d++) {
		if(fabs(total_deltaf[d]) > 1e-12) {
			printf(" step %d check_vel  momentum not conserved at (%d %d %d) (%le %le %le)!\n", nstep, i, j, k, total_deltaf[0], total_deltaf[1], total_deltaf[2]);
			break;    
		}
	}

}

void Write_time(int nstep, char *work_dir)
{
	extern int n_proc;
	char filename[160];
	FILE *stream;
  clock_t t = clock();

	sprintf(filename, "%s/data/run_time.dat", work_dir); 
	stream = fopen(filename, "a");
	fprintf(stream, "%d %d %le\n", n_proc, nstep, ((double)t)/CLOCKS_PER_SEC);
	fclose(stream);
}


void WallStress(int step, struct vector wallVel, Float ***velcs_df, char *work_dir, int
		outputInterval)
{
	// Wall stress should be calculated after fluid collision before fluid propagation.
	// wallVel: the velocity of the top moving wall

	if(step % outputInterval == 0)
	{
		extern int max_x, max_y, max_z;
		double wallArea = max_x*max_z;  //
		double sigmaYX = 0.;
		double viscosity = 0.0;
		double shearRate = 2.0 * wallVel.x / max_y;
		double strain = shearRate * step;
		struct vector deltaP={0., 0., 0.};
		char fileName[200];
		FILE *stream;
		sprintf(fileName, "%s/data/wallStress.dat", work_dir);
		stream = fopen(fileName, "a");

		for(int z=1; z <= max_z; z++) {  // loop fluid nodes on the top wall
			for(int x=1; x <= max_x; x++) {

				int xy = x * (max_y+2) + max_y;

				double fluidDensity = Rho_Fl;
//				for(int n=0; n < Num_Dir; n++) {
//					fluidDensity += velcs_df[xy][n][z];
//				}
				for(int n=0; n < Num_Dir_X; n++) {  // the components hitting the top wall
					int compon = y_p[n];
					double weighting = fac[compon] / 36.0;
					double movingBoundary = weighting * fluidDensity * (wallVel.x*c_x[compon] +
							wallVel.y*c_y[compon] + wallVel.z*c_z[compon]) / CS2;
					double density = 2.0 * (velcs_df[xy][compon][z] - movingBoundary);  // population density
					deltaP.x += density * c_x[compon];
					deltaP.y += density * c_y[compon];
					deltaP.z += density * c_z[compon];
				}
			}
		}
		// 1. Notice: the following calculation is correct only for dt=dx=1.
		// 2. The minus sign attribute to the normal vector of the top moving wall.
		deltaP.x /= (-wallArea);
		deltaP.y /= (-wallArea);
		deltaP.z /= (-wallArea);
		sigmaYX = deltaP.x;
		viscosity = sigmaYX / shearRate;
		fprintf(stream, "%lf %.4le %.4le %.4le %.4le\n",strain, viscosity, sigmaYX, 
        deltaP.y, deltaP.z); // (strain, viscosity, yx, yy, yz)
		fclose(stream);
	}
}

void WriteParticleStress(int step, int outputInterval, struct sphere_param *sphere_pm,
		 struct monomer *monomers, char *work_dir)
{
	if (step % outputInterval == 0)
	{
		extern int max_x, max_y, max_z;
		double boxSize = max_x*max_y*max_z;
		char fileName[200];
		FILE *stream;
		int beadNum = sphere_pm->num_beads;
		//double strain = shearRate * step;
		double stress_elas[3][3], stress_bend[3][3], stress_vol[3][3], stress_areaG[3][3], 
           stress_areaL[3][3], stress_wall[3][3], stress_inter_v1[3][3], 
           stress_inter_v2[3][3], stressElastic[3][3], stressAll[3][3];

		for(int i=0; i < DIMS; i++) {
			for(int j=0; j < DIMS; j++) {
        stress_elas[i][j]=0.;  stress_bend[i][j]=0.;  stress_vol[i][j]=0.;
        stress_areaG[i][j]=0.;  stress_areaL[i][j]=0.;
        stress_wall[i][j]=0.;
        stress_inter_v1[i][j]=0.;  stress_inter_v2[i][j]=0.;
        stressElastic[i][j]=0.;  stressAll[i][j]=0.;
      }
    }
		for(int n=0; n < beadNum; n++) {
			for(int i=0; i < DIMS; i++) {
				for(int j=0; j < DIMS; j++) {
          stress_elas[i][j]   += monomers[n].stress_elas[i][j];
          stress_bend[i][j]   += monomers[n].stress_bend[i][j];
          stress_vol[i][j]    += monomers[n].stress_vol[i][j];
          stress_areaG[i][j]  += monomers[n].stress_areaG[i][j];
          stress_areaL[i][j]  += monomers[n].stress_areaL[i][j];    
          stress_wall[i][j]   += monomers[n].stress_wall[i][j];
          stress_inter_v1[i][j] += monomers[n].stress_int_v1[i][j];
          stress_inter_v2[i][j] += monomers[n].stress_int_v2[i][j];

          stressElastic[i][j] += monomers[n].pos[i] *
                                 (monomers[n].force[j]-monomers[n].force_inter[j]);
          stressAll[i][j] += monomers[n].pos[i]*monomers[n].force[j];
				}
			}
		}
		for(int i=0; i < DIMS; i++) {
			for(int j=0; j < DIMS; j++) {
        stress_elas[i][j]   /= (-boxSize);
        stress_bend[i][j]   /= (-boxSize);
        stress_vol[i][j]    /= (-boxSize);
        stress_areaG[i][j]  /= (-boxSize);
        stress_areaL[i][j]  /= (-boxSize);
        stress_wall[i][j]   /= (-boxSize);
        stress_inter_v1[i][j] /= (-boxSize);
        stress_inter_v2[i][j] /= (-boxSize);
        
        stressElastic[i][j] /= (-boxSize);
        stressAll[i][j]     /= (-boxSize);
      }
    }
    double stress_v1[3][3], stress_v2[3][3];
    for(int i=0; i<3; i++) {
      for(int j=0; j<3; j++) {
        stress_v1[i][j] = stressElastic[i][j] + stress_inter_v1[i][j];
        stress_v2[i][j] = stressElastic[i][j] + stress_inter_v2[i][j];
      }
    }

 		sprintf(fileName,"%s/data/pStressElement.dat", work_dir);
		stream = fopen(fileName, "a");
		fprintf(stream, "%6d  ", step);
		fprintf(stream, "%+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le\n",
    stress_elas[1][0], stress_bend[1][0], stress_vol[1][0], stress_areaG[1][0], 
    stress_areaL[1][0], stress_wall[1][0], stress_inter_v1[1][0], stress_inter_v2[1][0]);
		fclose(stream);

    sprintf(fileName, "%s/data/pStressV1.dat", work_dir);
    stream = fopen(fileName, "a");
    fprintf(stream, "%6d  ", step);
    fprintf(stream, "%+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le\n",
        stress_v1[0][0],stress_v1[0][1],stress_v1[0][2],stress_v1[1][0],stress_v1[1][1],
        stress_v1[1][2],stress_v1[2][0],stress_v1[2][1],stress_v1[2][2]);
    fclose(stream);

    sprintf(fileName, "%s/data/pStressV2.dat", work_dir);
    stream = fopen(fileName, "a");
    fprintf(stream, "%6d  ", step);
    fprintf(stream, "%+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le\n",
        stress_v2[0][0],stress_v2[0][1],stress_v2[0][2],stress_v2[1][0],stress_v2[1][1],
        stress_v2[1][2],stress_v2[2][0],stress_v2[2][1],stress_v2[2][2]);
    fclose(stream);

    sprintf(fileName, "%s/data/pStressAll.dat", work_dir);
    stream = fopen(fileName, "a");
    fprintf(stream, "%6d  ", step);
    fprintf(stream, "%+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le\n",
        stressAll[0][0],stressAll[0][1],stressAll[0][2],stressAll[1][0],stressAll[1][1],
        stressAll[1][2],stressAll[2][0],stressAll[2][1],stressAll[2][2]);
    fclose(stream);

    sprintf(fileName, "%s/data/pStressElas.dat", work_dir);
    stream = fopen(fileName, "a");
    fprintf(stream, "%6d  ", step);
    fprintf(stream, "%+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le    %+.4le\n",
        stressElastic[0][0],stressElastic[0][1],stressElastic[0][2],stressElastic[1][0],stressElastic[1][1],
        stressElastic[1][2],stressElastic[2][0],stressElastic[2][1],stressElastic[2][2]);
    fclose(stream);

	}
}

void Write_Monomer(struct sphere *spheres, struct monomer *monomers, struct sphere_param 
		*sphere_pm, int n_step, int addsteps, char *work_dir)
{
	extern int max_x, max_y, max_z;
	static int buff_nstep[WRITE_BUFFER];
	static double buff_avgdisp2[WRITE_BUFFER], buff_avgdr2[WRITE_BUFFER], buff_avgvel[WRITE_BUFFER];
	static double buff_monvel[WRITE_BUFFER][500], buff_tempscale[WRITE_BUFFER];
	char filename[160];
	int j;
	int n, n2, d;
	int maxsize[DIMS];
	int num_write, num_step;
	int tot_time, add_time;
	double avg_disp2, avg_vel, avg_E, avg_dr2;
	int num_beads = sphere_pm->num_beads;
	FILE *stream;
	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;
	tot_time = (n_step+addsteps)*sphere_pm->MD_steps*sphere_pm->dt;
	add_time = addsteps*sphere_pm->MD_steps*sphere_pm->dt;

	if(n_step == 0) 
	{
		for(n=0; n<num_beads; n++) {
			for(d=0; d<DIMS; d++) {
				monomers[n].pos0[d]    = monomers[n].pos_pbc[d];
				monomers[n].pos_old[d] = monomers[n].pos_pbc[d];
			}
		}
		for(j=0; j<WRITE_BUFFER; j++) {
			buff_nstep[j]=0;
			buff_avgdisp2[j]=0.0;
			buff_avgdr2[j]=0.0;
			buff_avgvel[j]=0.0;
			buff_tempscale[j]=0.0;
		}
	}

//	if(sphere_pm->numsteps > 1001) 
//		//num_write = sphere_pm->relax_time/(sphere_pm->MD_steps*sphere_pm->dt);
//    num_write=5000;
//	else if(sphere_pm->num_beads == 1) 
//		num_write = 10;
//	else 
//		num_write = 1;
//
//	if((n_step+addsteps) % num_write == 0) 
//	{
//		/* Write out monomer position and velocity */
//		sprintf(filename, "%s/data/monpos.dat", work_dir);
//		stream=fopen(filename, "a");
//
//		if(num_beads < 500)
//		{ 
//			for(n=0; n<num_beads; n++)
//			{  
//				fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\t %12.10le %12.10le %12.10le\n", 	
//        tot_time,monomers[n].pos[0],monomers[n].pos[1],monomers[n].pos[2],
//        monomers[n].vel[0],monomers[n].vel[1],monomers[n].vel[2],
//        monomers[n].fluid_vel[0],monomers[n].fluid_vel[1],monomers[n].fluid_vel[2]);
//			}
//		}
//		else
//		{
//		  for(n=0; n<num_beads; n++) 
//				fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\n",
//						  tot_time,monomers[n].vel[0],monomers[n].vel[1],monomers[n].vel[2],
//              monomers[n].fluid_vel[0],monomers[n].fluid_vel[1],monomers[n].fluid_vel[2]);			
//		}  
//		fclose(stream);
//		/* write out bonding config */
//		if(WRITE_BOND)  // By default (head.h), WRITE_BOND==0!
//		{
//			sprintf(filename, "%s/data/bond.dat", work_dir);
//			stream=fopen(filename, "a");
//			for(n=0; n < num_beads ; n++) {  
//				for(j=1 ; j <= monomers[n].blist[0][0] ; j++) {
//					n2 = monomers[n].blist[j][0];
//					fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\n",
//							    tot_time, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2],
//							    monomers[n2].pos[0], monomers[n2].pos[1], monomers[n2].pos[2]);
//				}
//			}
//			fclose(stream);
//		}
//	}
//	/* write out monomer velocities */
//	if(num_beads < 5) {
//		for(n=0; n<num_beads; n++) {
//			sprintf(filename, "%s/data/monvel_%d.dat", work_dir, n);
//			stream=fopen(filename, "a");
//			fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le\n", tot_time, 
//              monomers[n].vel[0], monomers[n].vel[1], monomers[n].vel[2], 
//              monomers[n].fricforce[0]);
//			fclose(stream);
//		}
//  }
//20161129
	num_step = n_step / sphere_pm->write_time;
	if(n_step % sphere_pm->write_time == 0) 
  {
		avg_disp2 = 0.0;
		avg_vel = 0.0;
		avg_dr2 = 0.0;
    double avg_dr=0.; // 20161129

		for(n=0; n < num_beads; n++) 
    {
      double temp_disp=0.;
      double temp_dr=0.;
      double temp_vel=0.;
      monomers[n].dr2=0.;
			for(d=0; d < DIMS; d++) 
      {
        temp_disp = monomers[n].pos_pbc[d] - monomers[n].pos0[d];
				avg_disp2 += temp_disp * temp_disp;
        temp_dr =   monomers[n].pos_pbc[d] - monomers[n].pos_old[d];
        monomers[n].dr2 += temp_dr * temp_dr;
        temp_vel += monomers[n].vel[d] * monomers[n].vel[d];
				monomers[n].pos_old[d] = monomers[n].pos_pbc[d];
			}
      monomers[n].dr = sqrt(monomers[n].dr2);
      avg_dr += monomers[n].dr;
			avg_dr2 += monomers[n].dr2;
      avg_vel += sqrt(temp_vel); 
		}
		avg_disp2 /= num_beads;
    avg_dr /= num_beads;
		avg_dr2 /= num_beads;
    avg_vel /= num_beads;

		buff_nstep[num_step%WRITE_BUFFER] = tot_time;
		buff_avgdisp2[num_step%WRITE_BUFFER] = avg_disp2;
		buff_avgdr2[num_step%WRITE_BUFFER] = avg_dr; // avg_dr2; 20161129
		buff_avgvel[num_step%WRITE_BUFFER]=avg_vel;       //*Rho_Fl*sphere_pm->monmass;  
		buff_tempscale[num_step%WRITE_BUFFER] = sphere_pm->tempscale;  

		if((num_step+1)%WRITE_BUFFER == 0) {
			sprintf(filename, "%s/data/avg_disp2.dat", work_dir);
			stream = fopen(filename, "a");
			for(j=0; j<WRITE_BUFFER; j++)
				fprintf(stream, "%d %le %le %le\n", buff_nstep[j], buff_avgdisp2[j], buff_avgdr2[j],
                buff_avgvel[j]);
			fclose(stream);
		}
	} 
}

void WriteForce(int step, struct sphere_param *sphere_pm, struct monomer *mono, 
     char *work_dir)
{
  double elastic=0., bending=0., vol=0., areaG=0., areaL=0., total=0.;
  int numBead = sphere_pm->num_beads;
  for(int n=0; n < numBead; n++) {    
    elastic += sqrt(iproduct(mono[n].force_spring,mono[n].force_spring));
    bending += sqrt(iproduct(mono[n].force_bend,mono[n].force_bend));
    vol += sqrt(iproduct(mono[n].force_vol,mono[n].force_vol));
    areaG += sqrt(iproduct(mono[n].force_areaG,mono[n].force_areaG));
    areaL += sqrt(iproduct(mono[n].force_areaL,mono[n].force_areaL));
    total += sqrt(iproduct(mono[n].force,mono[n].force));
  }
  elastic/=numBead;
  bending/=numBead;
  vol/=numBead;
  areaG/=numBead;
  areaL/=numBead;
  total/=numBead;

  char filename[200];
  FILE *stream;
  sprintf(filename, "%s/data/vertexForce.dat", work_dir);
  stream = fopen(filename, "a"); 
  fprintf(stream, "%6d  ",step);
  fprintf(stream, "%+.4le    %+.4le    %+.4le    %+.4le   %+.4le    %+.4le\n",
  elastic, bending, vol, areaG, areaL, total);
  fclose(stream);
}

void WriteTemplate(struct sphere_param *sphere_pm, struct monomer *mono, char *work_dir)
{
  int numBead = sphere_pm->N_per_sphere[0];
  double com[DIMS]={0.};
  char filename[200];
  FILE *stream;
  sprintf(filename,"%s/init/n2_biconcave.dat", work_dir);
  double position[numBead][DIMS];

  for(int d=0; d < DIMS; d++) {
    for(int i=0; i < numBead; i++)
      com[d] += mono[i].pos_pbc[d];
    com[d] /= numBead;
  }
  for(int d=0; d < DIMS; d++) {
    for(int i=0; i < numBead; i++) {
      position[i][d] = mono[i].pos_pbc[d];
      position[i][d] -= com[d];
    }
  }

  stream = fopen(filename, "w");
  for(int n=0; n < numBead; n++)
    fprintf(stream, "%le %le %le\n", position[n][0], position[n][1], position[n][2]);
    //fprintf(stream, "%le %le %le\n", mono[n].pos_pbc[0], mono[n].pos_pbc[1], mono[n].pos_pbc[2]);
  for(int n=0; n < numBead; n++)
  {
    fprintf(stream, "%d %d\n", n, mono[n].blist[0][0]);
    for(int i=1; i <= mono[n].blist[0][0]; i++)
      fprintf(stream, "%d %d %d\n", mono[n].blist[i][0], mono[n].blist[i][1], 
             mono[n].blist[i][2]);
  }
  fclose(stream);
}

void BubbleSort(double *arr, int len)
{
  for(int i = 0; i < len - 1; i++)
    for(int j = 0; j < len - 1 - i; j++)
      if(arr[j] > arr[j + 1])
      {
        double temp = arr[j];
        arr[j] = arr[j + 1];
        arr[j + 1] = temp;
      }
}

void ParticleDim(struct sphere_param *sphere_pm, struct monomer *mono)
{
  int numBead = sphere_pm->num_beads;
  double dim[3][numBead];

  for(int d=0; d < DIMS; d++)
    for(int n=0; n < numBead; n++)
      dim[d][n] = mono[n].pos[d];
  for(int d=0; d < DIMS; d++)
    BubbleSort(dim[d], numBead);

  printf("particle dimensions (x,y,z)=(%f, %f, %f)\n",dim[0][numBead-1]-dim[0][0],
      dim[1][numBead-1]-dim[1][0], dim[2][numBead-1]-dim[2][0]);
}

// Functions below are not used.
//
/////////////////////////////////////////////////////////////////////////////////////////

void Discretize (int flow_flag, int window, int mark_interval, int oscillation_period, 
		int time_step, int **node_map, struct sphere_param *sphere_pm, 
		struct monomer *monomers, char *work_dir)
{
	static int index;
	static int temp_window;
	int half_window;
	int delta_interval;
	int print_interval;
	half_window    = window / 2;
	delta_interval = oscillation_period / 8;

	if (sphere_pm->num_beads != 0)
	{  
		if (flow_flag ==2)                            //  oscillatory shear
		{  

			print_interval = (oscillation_period / 8) - half_window + 1 + index*delta_interval;

			if (time_step > 0 && (time_step % print_interval== 0))
			{ 
				temp_window = window;
				index++;
			}

			if (temp_window > 0)
			{
				MarkParticle(time_step, node_map, sphere_pm, monomers, work_dir);
				temp_window--;
			}

		}
		if (flow_flag == 1)                           // unidirectional shear
		{ 
			if(time_step % mark_interval == 0)
			{
				temp_window = window ;
			}   

			if(temp_window > 0)
			{
				MarkParticle(time_step, node_map, sphere_pm, monomers, work_dir);
				temp_window--;
			}
		}  
	}
}

void Write_Particle(struct object *sphere, struct sphere_param *sphere_pm, 
		int n_step, char *work_dir)
{
	extern int max_x, max_y, max_z;
	extern int num_sph;

	static double x0[MAX_COLL], y0[MAX_COLL], z0[MAX_COLL];
	static int buff_nstep[WRITE_BUFFER];
	static double buff_avgdisp2[WRITE_BUFFER], buff_avgdr2[WRITE_BUFFER], buff_avgvel[WRITE_BUFFER];
	static double xp[MAX_COLL][WRITE_BUFFER], yp[MAX_COLL][WRITE_BUFFER], zp[MAX_COLL][WRITE_BUFFER];
	static double ux[MAX_COLL][WRITE_BUFFER], uy[MAX_COLL][WRITE_BUFFER], uz[MAX_COLL][WRITE_BUFFER];
	static double wx[MAX_COLL][WRITE_BUFFER], wy[MAX_COLL][WRITE_BUFFER], wz[MAX_COLL][WRITE_BUFFER];

	int i,j;
	int num_step;
	double dx, dy, dz;
	char filename[160];
	double avg_disp2;
	double avg_dr2;
	double avg_vel;
	FILE *stream;

	if(n_step == 0) 
	{
		for(i=0; i<num_sph; i++) 
		{
			x0[i]=sphere[i].r.x;
			y0[i]=sphere[i].r.y;
			z0[i]=sphere[i].r.z;
			sphere[i].rold.x=sphere[i].r.x;
			sphere[i].rold.y=sphere[i].r.y;
			sphere[i].rold.z=sphere[i].r.z;
			sphere[i].dr2 = 0.0;
		}

		for(j=0; j<WRITE_BUFFER; j++) 
		{
			buff_nstep[j]=0;
			buff_avgdisp2[j]=0.0;
			buff_avgvel[j]=0.0;

			for(i=0; i<MAX_COLL; i++) 
			{
				xp[i][j]=0.0; yp[i][j]=0.0; zp[i][j]=0.0;
				ux[i][j]=0.0; uy[i][j]=0.0; uz[i][j]=0.0;
				wx[i][j]=0.0; wy[i][j]=0.0; wz[i][j]=0.0;
			}
		}
	}

	num_step = n_step / sphere_pm->write_time;
	if(n_step % sphere_pm->write_time == 0) 
	{
		/* average mean squared displacement */
		for(i=0; i<num_sph; i++) 
		{
			dx = n_image(sphere[i].r.x - x0[i], max_x);
			dy = n_image(sphere[i].r.y - y0[i], max_y);
			dz = n_image(sphere[i].r.z - z0[i], max_z);
			sphere[i].disp2 = dx*dx+dy*dy+dz*dz;

			dx = n_image(sphere[i].r.x - sphere[i].rold.x, max_x);
			dy = n_image(sphere[i].r.y - sphere[i].rold.y, max_y);
			dz = n_image(sphere[i].r.z - sphere[i].rold.z, max_z);
			sphere[i].dr2 += dx*dx+dy*dy+dz*dz;
		}

		avg_disp2 = 0.0;
		avg_dr2 = 0.0;
		avg_vel = 0.0;
		for(i=0; i<num_sph; i++) 
		{
			avg_disp2 += sphere[i].disp2;
			avg_dr2 += sphere[i].dr2;
			avg_vel += sphere[i].u.x * sphere[i].u.x + sphere[i].u.y*sphere[i].u.y + sphere[i].u.z*sphere[i].u.z;
		}
		avg_disp2 /= num_sph;
		avg_dr2 /= num_sph;
		avg_vel /= num_sph;

		/* Store the data into the write buffer */
		buff_nstep[n_step%WRITE_BUFFER]=n_step;
		buff_avgdisp2[n_step%WRITE_BUFFER]=avg_disp2;
		buff_avgdr2[n_step%WRITE_BUFFER]=avg_dr2;
		buff_avgvel[n_step%WRITE_BUFFER]=avg_vel;
		for(i=0; i<num_sph; i++) 
		{
			xp[i][n_step%WRITE_BUFFER]=sphere[i].r.x;
			yp[i][n_step%WRITE_BUFFER]=sphere[i].r.y;
			zp[i][n_step%WRITE_BUFFER]=sphere[i].r.z;
			ux[i][n_step%WRITE_BUFFER]=sphere[i].u.x;
			uy[i][n_step%WRITE_BUFFER]=sphere[i].u.y;
			uz[i][n_step%WRITE_BUFFER]=sphere[i].u.z;
			wx[i][n_step%WRITE_BUFFER]=sphere[i].w.x;
			wy[i][n_step%WRITE_BUFFER]=sphere[i].w.y;
			wz[i][n_step%WRITE_BUFFER]=sphere[i].w.z;
		}

		for(i=0; i<num_sph; i++) 
		{
			sphere[i].rold.x=sphere[i].r.x;
			sphere[i].rold.y=sphere[i].r.y;
			sphere[i].rold.z=sphere[i].r.z;
		}

		if((num_step+1)%WRITE_BUFFER == 0) 
		{
			sprintf(filename, "%s/data/avg_disp2.dat", work_dir);
			stream = fopen(filename, "a");
			for(j=0; j<WRITE_BUFFER; j++)
				fprintf(stream, "%d %le %le %le\n", buff_nstep[j], buff_avgdisp2[j], buff_avgdr2[j], buff_avgvel[j]);
			fclose(stream);

			for(i=0; i<10; i++) 
			{
				sprintf(filename, "%s/data/coll%d.dat", work_dir, i);
				stream = fopen(filename, "a");
				for(j=0; j<WRITE_BUFFER; j++)
					fprintf(stream, "%d %le %le %le %le %le %le\n", buff_nstep[j], xp[i][j], yp[i][j], zp[i][j], ux[i][j], uy[i][j], uz[i][j]);
				fclose(stream);
			}
		}
	}
}

void Write_RDF(struct sphere *spheres, struct monomer *mon, 
		struct sphere_param *sphere_pm, int n_step, int addsteps, int num_step, 
		char *work_dir)
{
	extern int max_x, max_y, max_z;
	static int num_relax=0, num_write=0;
	static double gofy[BINS], gcofy[BINS];
	static int buff_gofy[WRITE_BUFFER][BINS][NTYPES], buff_gcofy[WRITE_BUFFER][BINS][NTYPES];
	static double buff_stretchofy[WRITE_BUFFER][BINS][NTYPES];
	static double buff_Rxofy[WRITE_BUFFER][BINS][NTYPES], buff_Ryofy[WRITE_BUFFER][BINS][NTYPES], buff_Rzofy[WRITE_BUFFER][BINS][NTYPES];
	static int buff_nstep[WRITE_BUFFER];

	char filename[160];
	int i,j, k;
	int n, d;
	int maxsize[DIMS];
	int Ntype[NTYPES];
	int Nbead[NTYPES];
	int num_beads[NTYPES];
	int tot_time, add_time, relax_step;
	int type, start;
	double *yc;
	//double range=(double)min(max_x, min(max_y, max_z));
	double range=max_y-1.;
	double binsize=range/BINS;
	FILE *stream;

	maxsize[0]=max_x;
	maxsize[1]=max_y;
	maxsize[2]=max_z;
	for (i=0 ; i<NTYPES ; i++) {
		Ntype[i] = sphere_pm->Ntype[i];
		Nbead[i] = sphere_pm->N_per_sphere[i];
		num_beads[i] = Ntype[i]*Nbead[i];
	}
	tot_time = (n_step+addsteps)*sphere_pm->MD_steps*sphere_pm->dt;
	add_time = addsteps*sphere_pm->MD_steps*sphere_pm->dt;
	relax_step = sphere_pm->relax_time/(sphere_pm->MD_steps*sphere_pm->dt);

	if((yc = (double *)malloc(sphere_pm->Nsphere*sizeof(double)))==NULL) exit(1);

	/* Calculate the time dependent g(r,t) */
	if(n_step == 0) {
		for(j=0; j<WRITE_BUFFER; j++) {
			buff_nstep[j]=0;
			for(i=0; i<BINS; i++)
				for (k=0 ; k<NTYPES ; k++){
					buff_gofy[j][i][k]=0;
					buff_gcofy[j][i][k]=0;
					buff_stretchofy[j][i][k]=0.0;
					buff_Rxofy[j][i][k]=0.0;
					buff_Ryofy[j][i][k]=0.0;
					buff_Rzofy[j][i][k]=0.0;
				}
		}
	}

	if(relax_step > 10) {
		if((n_step+addsteps) % (relax_step/10) == 0) {
			for(i=0; i<sphere_pm->Nsphere; i++) {
				type = (i<Ntype[0] ? 0 : 1);
				start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0]+(i-Ntype[0])*Nbead[1]);
				yc[i] = 0.0;
				for(j=0; j<Nbead[type]; j++) {
					yc[i]+=mon[start+j].pos[1];
					if (mon[start+j].pos[1] < range)
						buff_gofy[num_write%WRITE_BUFFER][(int)(mon[start+j].pos[1]/binsize)][type]++;
					else
						fprintf(stderr, "monomer %d y position out of range %le!\n", start+j, mon[start+j].pos[1]);
				}
				yc[i] /= Nbead[type];

				if(yc[i] < range) {
					buff_gcofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type]++;
					buff_stretchofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(spheres[i].stretch[0]*spheres[i].stretch[0]+spheres[i].stretch[1]*spheres[i].stretch[1]);
					buff_Rxofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(spheres[i].Rx2);
					buff_Ryofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(spheres[i].Ry2);
					buff_Rzofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(spheres[i].Rz2);
				}
				else
					fprintf(stderr, "sphere %d y position out of range %le!\n", i, yc[i]);
			}

			buff_nstep[num_write%WRITE_BUFFER]=tot_time;

			if((num_write+1)%WRITE_BUFFER == 0)
				for (k=0 ; k<NTYPES ; k++) {
					sprintf(filename, "%s/data/gofyoft-%d.dat", work_dir, k);
					stream=fopen(filename, "a");
					for(j=0; j<WRITE_BUFFER; j++)
						for(i=0; i<BINS; i++)
							fprintf(stream, "%d %le %le\n", buff_nstep[j], (double)i*binsize, (double)buff_gofy[j][i][k]/num_beads[k]);
					fclose(stream);

					sprintf(filename, "%s/data/gcofyoft-%d.dat", work_dir, k);
					stream=fopen(filename, "a");
					for(j=0; j<WRITE_BUFFER; j++)
						for(i=0; i<BINS; i++)
							fprintf(stream, "%d %le %le\n", buff_nstep[j], (double)i*binsize, (double)buff_gcofy[j][i][k]/Ntype[k]);
					fclose(stream);

					sprintf(filename, "%s/data/stretchofyoft-%d.dat", work_dir, k);
					stream=fopen(filename, "a");
					for(j=0; j<WRITE_BUFFER; j++)
						for(i=0; i<BINS; i++) {
							if (buff_gcofy[j][i][k] != 0) {
								buff_stretchofy[j][i][k] /= (double)buff_gcofy[j][i][k];
								buff_Rxofy[j][i][k] /= (double)buff_gcofy[j][i][k];
								buff_Ryofy[j][i][k] /= (double)buff_gcofy[j][i][k];
								buff_Rzofy[j][i][k] /= (double)buff_gcofy[j][i][k];
							}
							fprintf(stream, "%d %le %le %le %le %le\n", buff_nstep[j], (double)i*binsize, buff_stretchofy[j][i][k], buff_Rxofy[j][i][k], buff_Ryofy[j][i][k], buff_Rzofy[j][i][k]);
						}
					fclose(stream);

					for(j=0; j<WRITE_BUFFER; j++)
						for(i=0; i<BINS; i++) {
							buff_gofy[j][i][k]=0;
							buff_gcofy[j][i][k]=0;
							buff_stretchofy[j][i][k]=0.0;
							buff_Rxofy[j][i][k]=0.0;
							buff_Ryofy[j][i][k]=0.0;
							buff_Rzofy[j][i][k]=0.0;
						}
				}
			num_write++;
		}
	}  

	/* Calculate the overall g(r) and gc(r) */
	if(n_step ==0) {
		num_relax = addsteps/relax_step;
		printf("num_relax = %d\n", num_relax);

		for(i=0; i<BINS; i++) {
			gofy[i]=0.0;
			gcofy[i]=0.0;
		}

		if(num_relax > 0){
			sprintf(filename, "%s/data/gofy.dat", work_dir);
			if((stream = fopen(filename, "r")) != NULL) {
				fscanf(stream, "%*s %*s %*s %*s %*s %*s\n");
				for(i=0; i<BINS; i++)
					fscanf(stream, "%*s %lf\n", &gofy[i]);
				gofy[i] *= (double)num_relax;
				fclose(stream);
			}
			sprintf(filename, "%s/data/gcofy.dat", work_dir);
			if((stream = fopen(filename, "r")) != NULL) {
				for(i=0; i<BINS; i++)
					fscanf(stream, "%*s %lf\n", &gcofy[i]);
				gcofy[i] *= (double)num_relax;
				fclose(stream);
			}
		}
	}

	if((n_step+addsteps) % relax_step == 0) {
		for(i=0; i<sphere_pm->Nsphere; i++) {
			type = (i<Ntype[0] ? 0 : 1);
			start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0]+(i-Ntype[0])*Nbead[1]);
			yc[i] = 0.0;

			for(j=0; j<Nbead[type]; j++) {
				yc[i]+=mon[start+j].pos[1];

				gofy[(int)(mon[start+j].pos[1]/binsize)]++;
			}

			gcofy[(int)(yc[i]/(binsize*Nbead[type]))]++;
		}

		num_relax++;
	}

	if(n_step == num_step-1 && n_step > relax_step) {
		sprintf(filename, "%s/data/gofy.dat", work_dir);
		stream=fopen(filename, "w");
		fprintf(stream, "# %d beads %d Relaxation times \n", sphere_pm->num_beads, num_relax-1);
		for(i=0; i<BINS; i++)
			fprintf(stream, "%le %le\n", i*binsize, gofy[i]/(sphere_pm->num_beads*num_relax));
		fclose(stream);
		sprintf(filename, "%s/data/gcofy.dat", work_dir);
		stream=fopen(filename, "w");
		for(i=0; i<BINS; i++)
			fprintf(stream, "%le %le\n", i*binsize, gcofy[i]/(sphere_pm->Nsphere*num_relax));
		fclose(stream);
	}
	free(yc);
}
/* give eigenvales of a 3x3 symmetric matrix. algorithm taken from Smith (1961) */
void sym33eigen_values(double mx[3][3], double lambda[3])
{
	double m, p, q, phi, temp;
	int i, j;

	/* trace / 3 */
	m = (mx[0][0] + mx[1][1] + mx[2][2])/3.;

	mx[0][0] -= m;
	mx[1][1] -= m;
	mx[2][2] -= m;

	p = 0.;
	/* sum the square of elements */
	for (i=0 ; i<3 ; i++)
		for (j=0 ; j<3 ; j++)
			p+= mx[i][j]*mx[i][j];
	p /= 6.;

	/* determinant */
	q = mx[0][0]*mx[1][1]*mx[2][2] + mx[0][1]*mx[1][2]*mx[2][0] + mx[0][2]*mx[1][0]*mx[2][1]
		- (mx[0][0]*mx[1][2]*mx[2][1] + mx[0][1]*mx[1][0]*mx[2][2] + mx[0][2]*mx[1][1]*mx[2][0]);
	q /= 2.;

	temp = p*p*p - q*q;
	if (temp < 0.) {
		printf("rounding %le to 0!\n", temp);
		temp = 0.;
	}

	phi = atan(sqrt(temp)/q)/3.;

	lambda[0] = m + 2.*sqrt(p)*cos(phi);
	lambda[1] = m - sqrt(p)*(cos(phi)+sqrt(3.)*sin(phi));
	lambda[2] = m - sqrt(p)*(cos(phi)-sqrt(3.)*sin(phi));
}

void eigen_vector(double mx[3][3], double lambda, double v[3])
{
	int i;
	double a, b, c, d, delta, p, q, temp;

	a = mx[0][0] - lambda;
	b = mx[0][1];
	c = mx[1][0];
	d = mx[1][1] - lambda;
	delta = a*d - b*c;
	p = mx[0][2];
	q = mx[1][2];

	v[0] = (b*q-d*p)/delta;
	v[1] = (c*p-a*q)/delta;
	v[2] = 1.;

	temp = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	for (i=0;i<3;i++) v[i] /= temp;
}

void MarkParticle(int time_step, int **node_map, struct sphere_param *sphere_pm, struct 
		monomer *monomers, char *work_dir)
{
	extern int max_x, max_y, max_z;
	extern int num_x;
	int size_z = max_z+2;
	int size_y = max_y+2;
	int size_x = num_x+2;
	char filename[200];
	char name[200];
	FILE *stream;
	int node_num, bead_num;
	int i,j,k,n;

	sprintf(filename, "%s/data/mark_particle_%d.dat", work_dir, time_step);
	stream = fopen(filename, "w");

	node_num=0;
	for(k=0; k < size_z; k++)
	{
		for(j=0; j < size_y; j++)
		{
			for(i=0; i < size_x; i++)
			{
				int ij = i * size_y + j;

				if (node_map[ij][k]==2)          
				{	
					fprintf(stream, "%d %d %d\n", i,j,k);
					//	  printf("timestep = %d; node = %d\n", time_step, node_map[ij][k]);
					node_num++;
				}
			} 		  
		}
	}
	fclose(stream);      
	/*      
					sprintf(filename, "%s/data/mono_%d.dat", work_dir, time_step);
					stream = fopen(filename, "w");

					bead_num=0;
					for(n=0; n < sphere_pm->num_beads; n++)
					{
					fprintf(stream, "%lf %lf %lf\n", monomers[n].pos_pbc[0], monomers[n].pos_pbc[1],
					monomers[n].pos_pbc[2]);
					bead_num++; 
					}
					fclose(stream);
	 */
	sprintf(name, "%s/data/total_nodes_inside_%d.dat", work_dir, time_step);
	stream = fopen(name, "w");
	fprintf(stream, "%d\n", node_num/*, bead_num*/);
	fclose(stream);
}

void Write_stat(Float ***velcs_df, int nstep, double kT, char *work_dir)
{
	extern int max_x, max_y, max_z;
	extern int num_x;
	extern Float std_xi[Num_Dir];

	int i,j,k,q,p;
	int stat_obs[101][101][101];
	long nxy;
	char filename[160];
	double prop[Num_Dir];  
	static double stat_mom[101][101][101][Num_Dir];
	static double stat_momsq[101][101][101][Num_Dir];
	double ave_M, ave_sq_M, var_M;

	FILE *stream;

	// Statistics of the moments 

	/*  
			if(nstep == 100) {
			for(i=1; i<=num_x; i++)
			for(j=0; j<max_y; j++)
			for(k=0; k<max_z; k++) {
			stat_obs[i][j][k] = 0;

			for(q=0; q<10; q++) {
			stat_mom[i][j][k][q]=0.0;
			stat_momsq[i][j][k][q]=0.0;
			}
			}
			}

			if(nstep >=100) {
			for(i=1; i<=num_x; i++) 
			for(j=0; j<max_y; j++) {
			nxy = i*max_y + j;

			for(k=0; k<max_z; k++) {
			for(q=0; q < 10; q++)
			prop[q]=0.0;

			for(q=0; q < Num_Dir; q++) {
			prop[0] += evector[0][q]*velcs_df[nxy][q][k];
			prop[1] += evector[1][q]*velcs_df[nxy][q][k];
			prop[2] += evector[2][q]*velcs_df[nxy][q][k];
			prop[3] += evector[3][q]*velcs_df[nxy][q][k];
			prop[4] += evector[4][q]*velcs_df[nxy][q][k];
			prop[5] += evector[5][q]*velcs_df[nxy][q][k];
			prop[6] += evector[6][q]*velcs_df[nxy][q][k];
			prop[7] += evector[7][q]*velcs_df[nxy][q][k];
			prop[8] += evector[8][q]*velcs_df[nxy][q][k];
			prop[9] += evector[9][q]*velcs_df[nxy][q][k];

			}

			for(q=0; q < 10; q++) {
			stat_mom[i][j][k][q] +=prop[q];
			stat_momsq[i][j][k][q] +=prop[q]*prop[q];
			}

			stat_obs[i][j][k]++;
			}
			}
			}

			if(nstep == 1500) {
			for(q=0; q < 10; q++) {
			sprintf(filename, "%s/data/moment_%ld.dat", work_dir, q);

			stream = fopen(filename, "w");
			for(i=1; i<=num_x; i++) 
			for(j=0; j<max_y; j++) 
			for(k=0; k<max_z; k++) {
			ave_M = stat_mom[i][j][k][q]/(stat_obs[i][j][k]);
			ave_sq_M = stat_momsq[i][j][k][q]/(stat_obs[i][j][k]);
			var_M = ave_sq_M - ave_M*ave_M;

			if(q > 3)
			fprintf(stream, "%d %d %d\t%le\t%le\t%le\t%le\n", i,j,k, ave_M, ave_sq_M, var_M, var_M/(std_xi[q]*std_xi[q]));
			if(q <= 3)
			fprintf(stream, "%d %d %d\t%le\t%le\t%le\n", i,j,k, ave_M, ave_sq_M, var_M);

			}

			fclose(stream);
			}
			}

	 */

	for(i=1; i<=num_x; i++)
		for(j=1; j<=max_y; j++)
			for(k=1; k<=max_z; k++) 
				for(q=0; q<Num_Dir; q++) {
					stat_mom[i][j][k][q]=0.0;
					stat_momsq[i][j][k][q]=0.0;
				}

	for(i=1; i<=num_x; i++) 
		for(j=1; j<=max_y; j++) {
			nxy = i*(max_y+2) + j;

			for(k=1; k<=max_z; k++) {
				for(q=0; q < Num_Dir; q++)
					prop[q]=0.0;

				for(p=0; p < Num_Dir; p++)
					for(q=0; q < Num_Dir; q++) 
						prop[p] += evector[p][q]*velcs_df[nxy][q][k];

				for(q=0; q < Num_Dir; q++) {
					stat_mom[i][j][k][q] =prop[q];
					stat_momsq[i][j][k][q] =prop[q]*prop[q];
				}
			}
		}

	for(q=0; q < Num_Dir; q++) {
		ave_M = 0.0;
		ave_sq_M = 0.0;
		for(i=1; i<=num_x; i++) 
			for(j=1; j<=max_y; j++) 
				for(k=1; k<=max_z; k++) {
					ave_M += stat_mom[i][j][k][q];
					ave_sq_M += stat_momsq[i][j][k][q];
				}

		ave_M /= (num_x*max_y*max_z);     	  
		ave_sq_M /= (num_x*max_y*max_z); 
		var_M = ave_sq_M - ave_M*ave_M;

		sprintf(filename, "%s/data/moment_%ld.dat", work_dir, q);
		stream = fopen(filename, "a");

		if(q > 3)
			fprintf(stream, "%d\t%le\t%le\t%le\n", nstep, ave_M, ave_sq_M, var_M/(std_xi[q]*std_xi[q]));
		if(q <= 3)
			fprintf(stream, "%d\t%le\t%le\t%le\n", nstep, ave_M, ave_sq_M/Rho_Fl, var_M/Rho_Fl);

		fclose(stream);
	}
}

