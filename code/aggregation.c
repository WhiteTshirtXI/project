#include "header.h"

void WriteBlist(struct sphere_param *sphere_pm, struct monomer *monomers, char *work_dir)
{
	int index1=0;
	int index2=0;
	int numBead = sphere_pm->num_beads;
	char filename[200];
	FILE *stream;
	sprintf(filename, "%s/init/n2Blist.dat", work_dir);

	for(int i=0; i < numBead; i++)
		monomers[i].Blist[0][0] = 0;

	for(int i=0; i < numBead; i++) {
		for(int bond=1; bond <= monomers[i].blist[0][0]; bond++)
		{
			int j = monomers[i].blist[bond][0];
			index1 = ++monomers[i].Blist[0][0];
			index2 = ++monomers[j].Blist[0][0];
			monomers[i].Blist[index1][0] = monomers[i].blist[bond][0];
			monomers[i].Blist[index1][1] = monomers[i].blist[bond][1];
			monomers[i].Blist[index1][2] = monomers[i].blist[bond][2];
			monomers[j].Blist[index2][0] = i;
			monomers[j].Blist[index2][1] = monomers[i].blist[bond][2];
			monomers[j].Blist[index2][2] = monomers[i].blist[bond][1];
		}
	}
	stream = fopen(filename, "w");  
	for(int i=0; i < sphere_pm->N_per_sphere[0]; i++) 
  {
		fprintf(stream, "%d %d\n", i, monomers[i].Blist[0][0]);

		for(int bond=1; bond <= monomers[i].Blist[0][0]; bond++) 
			fprintf(stream, "%d %d %d\n", monomers[i].Blist[bond][0], 
          monomers[i].Blist[bond][1], monomers[i].Blist[bond][2]);
	}
	fclose(stream);
} 

void AssignBlist(struct sphere_param *sphere_pm, struct monomer *monomers, char *work_dir)
{
	int Nbeads = sphere_pm->N_per_sphere[0];
	int Nspheres = sphere_pm->Ntype[0];
	int num_beads = sphere_pm->num_beads;
	char filename[200];
	FILE *stream;
  /*
	
		 if (sphere_pm->nlevel[0] == -1)       // (syfan 151125)
		 sprintf(filename, "%s/init/n1_rbc_Blist.init", work_dir);
		 else if (sphere_pm->nlevel[0] == -2)
		 sprintf(filename, "%s/init/n2_rbc_Blist.init", work_dir);
		 else if (sphere_pm->nlevel[0] == -3)
		 sprintf(filename, "%s/init/n3_rbc_Blist.init", work_dir);
		 else if (sphere_pm->nlevel[0] == -4)
		 sprintf(filename, "%s/init/n4_rbc_Blist.init", work_dir);
	*/ 
  if(sphere_pm->nlevel[0]==2)
    sprintf(filename, "%s/init/n2Blist.dat", work_dir);
  if(sphere_pm->nlevel[0]==3)
    sprintf(filename, "%s/init/n3Blist.dat", work_dir);
  if(sphere_pm->nlevel[0]==-1)
    sprintf(filename, "%s/init/n2Blist.dat", work_dir);  // 20170307 temporary
	stream = fopen(filename, "r");
	for(int i=0; i < Nbeads; i++) 
	{
		fscanf(stream, "%*s %d\n", &monomers[i].Blist[0][0]);
		for(int j=1; j <= monomers[i].Blist[0][0]; j++)
			fscanf(stream, "%d %d %d\n", &monomers[i].Blist[j][0], &monomers[i].Blist[j][1], 
					&monomers[i].Blist[j][2]);
	}
	fclose(stream);

	for(int i = 1; i < Nspheres; i++) 
	{
		int start = i * Nbeads;
		for(int j=0; j < Nbeads; j++) 
		{
			monomers[start + j].Blist[0][0] = monomers[j].Blist[0][0];
			for(int n=1; n <= monomers[j].Blist[0][0]; n++) {
				monomers[start + j].Blist[n][0] = start + (monomers[j].Blist[n][0]);
				monomers[start + j].Blist[n][1] = start + (monomers[j].Blist[n][1]);
				monomers[start + j].Blist[n][2] = start + (monomers[j].Blist[n][2]);
			}
		}
	}
}

void SetFace(struct sphere_param *sphere_pm, struct monomer *monomers, struct face *faces)
{
	int Ntype[NTYPES];
	int Nbeads[NTYPES];
	int Nfaces[NTYPES];
	int num_beads = sphere_pm->num_beads;
	int nface = 0;
	int numBond;
	int i, j, k;
	int a, b, c;
	int aa, bb, cc;
	unsigned check_face = 0;
	for (i = 0; i < NTYPES; i++) {
		Ntype[i] = sphere_pm->Ntype[i];
		Nbeads[i] = sphere_pm->N_per_sphere[i];
		Nfaces[i] = sphere_pm->face_per_sphere[i];
	}

	for(i=0; i < num_beads; i++) 
  {
		numBond = monomers[i].Blist[0][0];
		for(j=1; j <= numBond; j++) 
		{
			a = i;
			b = monomers[i].Blist[j][0];
			check_face = 0;
			/* Check if the face is existing, and set the face list  */
			for(k = 0; k < nface; k++) 
      {
				aa = faces[k].vertices[0];
				bb = faces[k].vertices[1];
				cc = faces[k].vertices[2];
				/* the three vertices of face mon[i].flist[j] is (i, mon[i].Blist[j][0], mon[i].Blist[j][1]) */
				if (bb == a && cc == b) {
					check_face = 1;
					monomers[i].faceList[j] = k;
					break;
				}
				else if (cc == a && aa == b) {
					check_face = 1;
					monomers[i].faceList[j] = k;
					break;
				}
			}
			/* Added a new face */
			if(check_face == 0) 
      {
				if(i < Ntype[0] * Nbeads[0])
					faces[nface].sphere_id = i / Nbeads[0];
				else
					faces[nface].sphere_id = Ntype[0] + (i - Ntype[0] * Nbeads[0]) / Nbeads[1];

				faces[nface].vertices[0] = i;
				faces[nface].vertices[1] = monomers[i].Blist[j][0];
				faces[nface].vertices[2] = monomers[i].Blist[j][1];
				monomers[i].faceList[j] = nface;
				nface++;
			}
		}
	}
	printf("%d faces added. %d faces / sphere0  %d faces / sphere1, total faces %d\n", 
	       nface, Nfaces[0], Nfaces[1], Nfaces[0] * Ntype[0] + Nfaces[1] * Ntype[1]);
	if(nface != Nfaces[0] * Ntype[0] + Nfaces[1] * Ntype[1])
	  fatal_err("Number of faces does not match", -1);
}

int* Sort(double arr[], int len, int *index) 
{
  for(int i=0; i < len; i++)
    index[i] = i+1;

	for(int i=0; i < len - 1; i++)
		for(int j=0; j < len - 1 - i; j++)
			if(arr[j] > arr[j + 1]) 
      {
/*
				temp = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = temp;
*/
        int temp_index = index[j];
        index[j] = index[j+1];
        index[j+1] = temp_index;
			}
  return index;
}

void sort_label(double arr[], int *index) 
{
  int len = sizeof(index)/sizeof(int);
	for(int i=0; i < len-1; i++)
		for(int j=0; j < len-1-i; j++)
			if(arr[j] > arr[j+1]) 
      {
        int temp_index = index[j];
        index[j] = index[j+1];
        index[j+1] = temp_index;
			}
}

//void AggreForce(struct sphere_param *sphere_pm, struct monomer *monomers)
//{  
//  int totalBead = sphere_pm->num_beads;
//  double sigma = sphere_pm->range_LJ;
//  double cutoff = 1.122 * sigma;
//  double epsWCA = sphere_pm->eps_overlap * sphere_pm->kT; 
//  double truncation = 1.0 * sigma;
//
//  for(int n1=0; n1 < totalBead; n1++) 
//  { 
//    int numNeighbor;
//    int nearest;
//    double dis_nearest;
//    int numBond, secBond, thirdBond, second, third;
//    double v12[DIMS];
//    double faceNormal[DIMS];
//    double faceCOM[DIMS];
//    double temp1[DIMS];
//    double temp2[DIMS];
//    double temp3[DIMS];
//    double dis1_inverse=0.;
//    double dis2_inverse=0.;
//    double dis3_inverse=0.;
//    double dis_inverse=0.;
//    double dis_beadFace=0.;
//    double aterm=0.;  
//    double force[DIMS]={0.0};
//    extern int max_x, max_y, max_z;
//    extern int wall_flag;
//    int maxsize[DIMS];
//    numNeighbor = monomers[n1].list[0];
//    dis_nearest = 100.0;
//    maxsize[0]=max_x;
//    maxsize[1]=max_y;
//    maxsize[2]=max_z;
//    int onSameParticle=1;
//
////printf("numNeighber=%d\n", numNeighbor);
//
//    // 1. find the nearest neighbor, on the other particle, in n1's neighbor list.      
//    for(int index=1; index <= numNeighbor; index++) 
//    {             
//      double dis_temp=0.0;
//      int n2 = monomers[n1].list[index];
//      if(monomers[n1].sphere_id != monomers[n2].sphere_id) // on different particles
//      {
//        onSameParticle=0;
//
//        for(int j=0; j < DIMS; j++) {
//          v12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
//          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1)) 
//            v12[j] = n_image( v12[j], maxsize[j] );
//          dis_temp += v12[j]*v12[j];
//        }
//        dis_temp = sqrt(dis_temp);
//
//        if(dis_temp < dis_nearest) {
//          nearest = n2;
//          dis_nearest = dis_temp;
//        }
//      }
//    }// get the nearest neighbor's index and the distance bwtween it.       
//    // 2. Call nearest neighbor's Bond list to find n1's top 3 nearest neighbors, including n2.
//    
//    if(onSameParticle==0)
//    {
//      numBond = monomers[nearest].Blist[0][0];
//      double dis[numBond+1];
//      double dis_temp[numBond];
//      int rank[numBond];
//      double bond0[DIMS]; 
//      double bond[DIMS];
//      double norm2=0.;
//      dis[0] = dis_nearest;
//
////printf("numBond=%d\n", numBond);
////printf("dis=(%f %f %f)\n",dis[0],dis[1],dis[2]);
////printf("dis_temp=(%f, %f, %f)\n",dis_temp[0],dis_temp[1],dis_temp[2]);
////printf("rank=(%d %d %d)\n", rank[0],rank[1],rank[2]);
//
//      for(int i=1; i <= numBond; i++)
//      {
//        int index = monomers[nearest].Blist[i][0];
//
//        for(int j=0; j < DIMS; j++) {
//          v12[j] = monomers[n1].pos_pbc[j] - monomers[index].pos_pbc[j];  // relative position vector (n1-index)
//          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
//            v12[j] = n_image(v12[j], maxsize[j]);
//          dis[i] += v12[j]*v12[j];
//        }
//        dis[i] = sqrt(dis[i]);
//      }
//      for(int i=0; i < numBond; i++)
//        dis_temp[i] = dis[i+1];
//
//      Sort(dis_temp, numBond, rank);    
//      secBond = rank[0];   // the bond index for the second nearest neighbor of n1
//      thirdBond = rank[1]; // the bond index for the third nearest neighbor of n1
//      second = monomers[nearest].Blist[secBond][0];  // the monomer index of the 2nd nearest neighbor
//      third = monomers[nearest].Blist[thirdBond][0];
///*
//printf("numBond=%d\n", numBond);
//printf("dis[0]=%f dis[1]=%f\n", dis[0],dis[1]);
//printf("secBond=%d  thirdBond=%d  second=%d  third=%d\n", secBond, thirdBond, second, third);
//PAUSE
//*/
//      // 3. Calculate the face's normal vector, the relative dist., and the interation
//      if(monomers[nearest].Blist[secBond][1] == third || 
//         monomers[nearest].Blist[secBond][2] == third) // the three vertexes belong to a face
//      {
//        // 3.1 normal vector 
//        for(int d=0; d < DIMS; d++) {
//          bond0[d] = monomers[second].pos_pbc[d] - monomers[nearest].pos_pbc[d];
//          bond[d] = monomers[third].pos_pbc[d] - monomers[nearest].pos_pbc[d];
//        }
//        if(monomers[nearest].Blist[secBond][1] == third) {
//          product(bond, bond0, faceNormal); 
//        }
//        else {
//          product(bond0, bond, faceNormal);
//        }
//        norm2 = iproduct(faceNormal, faceNormal);
//        for(int d=0; d < DIMS; d++)
//          faceNormal[d] = faceNormal[d] / sqrt(norm2);
//
//        // 3.2 face COM
//        for(int d=0; d < DIMS; d++)
//          faceCOM[d] = (monomers[nearest].pos_pbc[d] + monomers[second].pos_pbc[d] + 
//                        monomers[third].pos_pbc[d]) / 3.;
//        // distances between faceCOM and 3 vertexes 
//        for(int d=0; d < DIMS; d++) {
//          temp1[d] = faceCOM[d] - monomers[nearest].pos_pbc[d];
//          temp2[d] = faceCOM[d] - monomers[second].pos_pbc[d];
//          temp3[d] = faceCOM[d] - monomers[third].pos_pbc[d];
//          dis1_inverse += temp1[d]*temp1[d];
//          dis2_inverse += temp2[d]*temp2[d];
//          dis3_inverse += temp3[d]*temp3[d];
//        }
//        dis1_inverse = 1./ sqrt(dis1_inverse);
//        dis2_inverse = 1./ sqrt(dis2_inverse);
//        dis3_inverse = 1./ sqrt(dis3_inverse);
//        dis_inverse = dis1_inverse + dis2_inverse + dis3_inverse;
//
//        // 3.2. the relative distance between n1 and the face
//        for(int d=0; d < DIMS; d++)
//          v12[d] = monomers[n1].pos_pbc[d] - monomers[nearest].pos_pbc[d];
//        dis_beadFace = iproduct(v12, faceNormal);
//
//        // 3.3 L.-J. interaction
//        if(dis_beadFace < cutoff)
//        {
//          if(dis_beadFace < truncation)
//          {
////printf("within truncation range\n");
//            aterm = (sigma / truncation);
//            aterm = aterm * aterm * aterm * aterm * aterm * aterm;
//            for(int d = 0; d < DIMS; d++)
//              force[d] = 24.0 * (epsWCA/truncation) * (2.0*(aterm*aterm) - aterm) * faceNormal[d];
//            for(int d=0; d < DIMS; d++) {
//              monomers[n1].force[d] += force[d];
//              monomers[n1].force_face[d] += force[d];
//
//              monomers[nearest].force[d] -= (force[d] * dis1_inverse / dis_inverse);
//              monomers[second].force[d] -= (force[d] * dis2_inverse / dis_inverse);
//              monomers[third].force[d] -= (force[d] * dis3_inverse / dis_inverse);
//              monomers[nearest].force_face[d] -= (force[d] * dis1_inverse / dis_inverse);
//              monomers[second].force_face[d] -= (force[d] * dis2_inverse / dis_inverse);
//              monomers[third].force_face[d] -= (force[d] * dis3_inverse / dis_inverse);
//            }
//          }
//          else { 
////printf("interact\n");
//            aterm = (sigma / dis_beadFace);
//            aterm = aterm * aterm * aterm * aterm * aterm * aterm;
//            for(int d = 0; d < DIMS; d++)
//              force[d] = 24.0 * (epsWCA/dis_beadFace) * (2.0*(aterm*aterm) - aterm) * faceNormal[d];
//            for(int d=0; d < DIMS; d++) {
//              monomers[n1].force[d] += force[d];
//              monomers[n1].force_face[d] += force[d];
//
//              monomers[nearest].force[d] -= (force[d] * dis1_inverse / dis_inverse);
//              monomers[second].force[d] -= (force[d] * dis2_inverse / dis_inverse);
//              monomers[third].force[d] -= (force[d] * dis3_inverse / dis_inverse);
//              monomers[nearest].force_face[d] -= (force[d] * dis1_inverse / dis_inverse);
//              monomers[second].force_face[d] -= (force[d] * dis2_inverse / dis_inverse);
//              monomers[third].force_face[d] -= (force[d] * dis3_inverse / dis_inverse);
//            }
//          }
//        }   
///*
//        else {
//          for(int d=0; d < DIMS; d++) {
//            force[d] = 0.;
//            monomers[n1].force[d] += force[d];
//            monomers[n1].force_face[d] += force[d];
//          }
//        }
//*/                
//      }
//      /* 
//      else
//      {
//          /*
//          method 1:
//          -distance: bead-to-bead
//          -direction: the average normal vector of faces around the bead
//          -object: n1  and beads in n2's bond list
//         
//          method 2:
//          -object: n1 and beads on the nearest face
//          -distance: the distance between n1 and the face
//          -direction: the normal vector of the face
//          */
//    //}
//    }            
//  }
//}
void AggreForce(struct sphere_param *sphere_pm, struct monomer *monomers)
{  
  int totalBead = sphere_pm->num_beads;
//  double sigma = sphere_pm->range_LJ;
  double sigma = sphere_pm->eq_LJ / pow(2., 1./6.);
  double cutoff = 1.122 * sigma;
  double epsWCA = sphere_pm->epsilon_overlap * sphere_pm->kT; 
  double truncation = 1.0 * sigma;

  for(int n1=0; n1 < totalBead; n1++) 
  { 
    int numNeighbor;
    int nearest;
    double dis_nearest;
    int numBond, secBond, thirdBond, second, third;
    double v12[DIMS];
    double faceNormal[DIMS];
    double faceCOM[DIMS];
    double temp1[DIMS];
    double temp2[DIMS];
    double temp3[DIMS];
    double dis1_inverse=0.;
    double dis2_inverse=0.;
    double dis3_inverse=0.;
    double dis_inverse=0.;
    double dis_beadFace=0.;
    double aterm=0.;  
    double force[DIMS]={0.0};
    extern int max_x, max_y, max_z;
    extern int wall_flag;
    int maxsize[DIMS];
    numNeighbor = monomers[n1].list[0];
    dis_nearest = 100.0;
    maxsize[0]=max_x;
    maxsize[1]=max_y;
    maxsize[2]=max_z;
    int onSameParticle=1;

//printf("numNeighber=%d\n", numNeighbor);

    // 1. find the nearest neighbor, on the other particle, in n1's neighbor list.      
    for(int index=1; index <= numNeighbor; index++) 
    {             
      double dis_temp=0.0;
      int n2 = monomers[n1].list[index];
      if(monomers[n1].sphere_id != monomers[n2].sphere_id) // on different particles
      {
        onSameParticle=0;

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1)) 
            v12[j] = n_image( v12[j], maxsize[j] );
          dis_temp += v12[j]*v12[j];
        }
        dis_temp = sqrt(dis_temp);

        if(dis_temp < dis_nearest) {
          nearest = n2;
          dis_nearest = dis_temp;
        }
      }
    }// get the nearest neighbor's index and the distance bwtween it.       
    // 2. Call nearest neighbor's Bond list to find n1's top 3 nearest neighbors, including n2.
    
    if(onSameParticle==0)
    {
      numBond = monomers[nearest].Blist[0][0];
      double dis[numBond+1];
      double dis_temp[numBond];
      int rank[numBond];
      double bond0[DIMS]; 
      double bond[DIMS];
      double norm2=0.;
      dis[0] = dis_nearest;

//printf("numBond=%d\n", numBond);
//printf("dis=(%f %f %f)\n",dis[0],dis[1],dis[2]);
//printf("dis_temp=(%f, %f, %f)\n",dis_temp[0],dis_temp[1],dis_temp[2]);
//printf("rank=(%d %d %d)\n", rank[0],rank[1],rank[2]);

      for(int i=1; i <= numBond; i++)
      {
        int index = monomers[nearest].Blist[i][0];

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[index].pos_pbc[j];  // relative position vector (n1-index)
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
            v12[j] = n_image(v12[j], maxsize[j]);
          dis[i] += v12[j]*v12[j];
        }
        dis[i] = sqrt(dis[i]);
      }
      for(int i=0; i < numBond; i++)
        dis_temp[i] = dis[i+1];

      Sort(dis_temp, numBond, rank);    
      secBond = rank[0];   // the bond index for the second nearest neighbor of n1
      thirdBond = rank[1]; // the bond index for the third nearest neighbor of n1
      second = monomers[nearest].Blist[secBond][0];  // the monomer index of the 2nd nearest neighbor
      third = monomers[nearest].Blist[thirdBond][0];
/*
printf("numBond=%d\n", numBond);
printf("dis[0]=%f dis[1]=%f\n", dis[0],dis[1]);
printf("secBond=%d  thirdBond=%d  second=%d  third=%d\n", secBond, thirdBond, second, third);
PAUSE
*/
      // 3. Calculate the face's normal vector, the relative dist., and the interation
      if(monomers[nearest].Blist[secBond][1] == third || 
         monomers[nearest].Blist[secBond][2] == third) // the three vertexes belong to a face
      {
        // 3.1 normal vector 
        for(int d=0; d < DIMS; d++) {
          bond0[d] = monomers[second].pos_pbc[d] - monomers[nearest].pos_pbc[d];
          bond[d] = monomers[third].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        }
        if(monomers[nearest].Blist[secBond][1] == third) {
          product(bond, bond0, faceNormal); 
        }
        else {
          product(bond0, bond, faceNormal);
        }
        norm2 = iproduct(faceNormal, faceNormal);
        for(int d=0; d < DIMS; d++)
          faceNormal[d] = faceNormal[d] / sqrt(norm2);

        // 3.2 face COM
        for(int d=0; d < DIMS; d++)
          faceCOM[d] = (monomers[nearest].pos_pbc[d] + monomers[second].pos_pbc[d] + 
                        monomers[third].pos_pbc[d]) / 3.;
        // distances between faceCOM and 3 vertexes 
        for(int d=0; d < DIMS; d++) {
          temp1[d] = faceCOM[d] - monomers[nearest].pos_pbc[d];
          temp2[d] = faceCOM[d] - monomers[second].pos_pbc[d];
          temp3[d] = faceCOM[d] - monomers[third].pos_pbc[d];
          dis1_inverse += temp1[d]*temp1[d];
          dis2_inverse += temp2[d]*temp2[d];
          dis3_inverse += temp3[d]*temp3[d];
        }
        dis1_inverse = 1./ sqrt(dis1_inverse);
        dis2_inverse = 1./ sqrt(dis2_inverse);
        dis3_inverse = 1./ sqrt(dis3_inverse);
        dis_inverse = dis1_inverse + dis2_inverse + dis3_inverse;

        // 3.2. the relative distance between n1 and the face
        for(int d=0; d < DIMS; d++)
          v12[d] = monomers[n1].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        dis_beadFace = iproduct(v12, faceNormal);

        // 3.3 L.-J. interaction
        if(dis_beadFace < cutoff)
        {
          if(dis_beadFace < truncation)
          {
//printf("within truncation range\n");
            aterm = (sigma / truncation);
            aterm = aterm * aterm * aterm * aterm * aterm * aterm;
            for(int d = 0; d < DIMS; d++)
              force[d] = 24.0 * (epsWCA/truncation) * (2.0*(aterm*aterm) - aterm) * faceNormal[d];
            for(int d=0; d < DIMS; d++) {
              monomers[n1].force[d] += force[d];
              monomers[n1].force_face[d] += force[d];

              monomers[nearest].force[d] -= (force[d] * dis1_inverse / dis_inverse);
              monomers[second].force[d] -= (force[d] * dis2_inverse / dis_inverse);
              monomers[third].force[d] -= (force[d] * dis3_inverse / dis_inverse);
              monomers[nearest].force_face[d] -= (force[d] * dis1_inverse / dis_inverse);
              monomers[second].force_face[d] -= (force[d] * dis2_inverse / dis_inverse);
              monomers[third].force_face[d] -= (force[d] * dis3_inverse / dis_inverse);
            }
          }
          else { 
//printf("interact\n");
            aterm = (sigma / dis_beadFace);
            aterm = aterm * aterm * aterm * aterm * aterm * aterm;
            for(int d = 0; d < DIMS; d++)
              force[d] = 24.0 * (epsWCA/dis_beadFace) * (2.0*(aterm*aterm) - aterm) * faceNormal[d];
            for(int d=0; d < DIMS; d++) {
              monomers[n1].force[d] += force[d];
              monomers[n1].force_face[d] += force[d];

              monomers[nearest].force[d] -= (force[d] * dis1_inverse / dis_inverse);
              monomers[second].force[d] -= (force[d] * dis2_inverse / dis_inverse);
              monomers[third].force[d] -= (force[d] * dis3_inverse / dis_inverse);
              monomers[nearest].force_face[d] -= (force[d] * dis1_inverse / dis_inverse);
              monomers[second].force_face[d] -= (force[d] * dis2_inverse / dis_inverse);
              monomers[third].force_face[d] -= (force[d] * dis3_inverse / dis_inverse);
            }
          }
        }   
/*
        else {
          for(int d=0; d < DIMS; d++) {
            force[d] = 0.;
            monomers[n1].force[d] += force[d];
            monomers[n1].force_face[d] += force[d];
          }
        }
*/                
      }
      /* 
      else
      {
          /*
          method 1:
          -distance: bead-to-bead
          -direction: the average normal vector of faces around the bead
          -object: n1  and beads in n2's bond list
         
          method 2:
          -object: n1 and beads on the nearest face
          -distance: the distance between n1 and the face
          -direction: the normal vector of the face
          */
    //}
    }            
  }
}

void AggreForce2(struct sphere_param *sphere_pm, struct monomer *monomers, char *work_dir)
{  
  int totalBead = sphere_pm->num_beads;
//  double sigma = sphere_pm->range_LJ;
//  double sigma = sphere_pm->eq_LJ / pow(2., 1./6.);
//  double cutoff = 1.122 * sigma;
  double cutoff = sphere_pm->eq_LJ;
  char filename[100];
  FILE *stream;
  sprintf(filename, "%s/data/dis_beadFace.dat", work_dir);
  stream = fopen(filename, "a");

  for(int n1=0; n1 < totalBead; n1++) 
  { 
    int numNeighbor;
    int nearest;
    double dis_nearest;
    int numBond, secBond, thirdBond, second, third;
    double v12[DIMS];
    double faceNormal[DIMS];
    double dis_beadFace=0.;
    extern int max_x, max_y, max_z;
    extern int wall_flag;
    int maxsize[DIMS];
    numNeighbor = monomers[n1].list[0];
    dis_nearest = 100.0;
    maxsize[0]=max_x;
    maxsize[1]=max_y;
    maxsize[2]=max_z;
    int onSameParticle=1;
//printf("numNeighber=%d\n", numNeighbor);

    // 1. find the nearest neighbor, on the other particle, in n1's neighbor list.      
    for(int index=1; index <= numNeighbor; index++) 
    {             
      double dis_temp=0.0;
      int n2 = monomers[n1].list[index];
      if(monomers[n1].sphere_id != monomers[n2].sphere_id) // on different particles
      {
        onSameParticle=0;

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1)) 
            v12[j] = n_image( v12[j], maxsize[j] );
          dis_temp += v12[j]*v12[j];
        }
        dis_temp = sqrt(dis_temp);

        if(dis_temp < dis_nearest) {
          nearest = n2;
          dis_nearest = dis_temp;
        }
      }
    }// get the nearest neighbor's index and the distance bwtween it.       

    // 2. Call nearest neighbor's Bond list to find n1's top 3 nearest neighbors, including n2.    
    if(onSameParticle==0)
    {
      numBond = monomers[nearest].Blist[0][0];
      double dis[numBond+1];
      double dis_temp[numBond];
      int rank[numBond];
      double bond0[DIMS]; 
      double bond[DIMS];
      double norm2=0.;
      dis[0] = dis_nearest;
//printf("numBond=%d\n", numBond);
//printf("dis=(%f %f %f)\n",dis[0],dis[1],dis[2]);
//printf("dis_temp=(%f, %f, %f)\n",dis_temp[0],dis_temp[1],dis_temp[2]);
//printf("rank=(%d %d %d)\n", rank[0],rank[1],rank[2]);

      for(int i=1; i <= numBond; i++)
      {
        int index = monomers[nearest].Blist[i][0];

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[index].pos_pbc[j];  // relative position vector (n1-index)
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
            v12[j] = n_image(v12[j], maxsize[j]);
          dis[i] += v12[j]*v12[j];
        }
        dis[i] = sqrt(dis[i]);
      }
      for(int i=0; i < numBond; i++)
        dis_temp[i] = dis[i+1];

      Sort(dis_temp, numBond, rank);    
      secBond = rank[0];   // the bond index for the second nearest neighbor of n1
      thirdBond = rank[1]; // the bond index for the third nearest neighbor of n1
      second = monomers[nearest].Blist[secBond][0];  // the monomer index of the 2nd nearest neighbor
      third = monomers[nearest].Blist[thirdBond][0];
/*
printf("numBond=%d\n", numBond);
printf("dis[0]=%f dis[1]=%f\n", dis[0],dis[1]);
printf("secBond=%d  thirdBond=%d  second=%d  third=%d\n", secBond, thirdBond, second, third);
PAUSE
*/
      // 3. Calculate the face's normal vector, the relative dist., and the interation
      if(monomers[nearest].Blist[secBond][1] == third || 
         monomers[nearest].Blist[secBond][2] == third) // the three vertexes belong to a face
      {
        // 3.1 normal vector 
        for(int d=0; d < DIMS; d++) {
          bond0[d] = monomers[second].pos_pbc[d] - monomers[nearest].pos_pbc[d];
          bond[d] = monomers[third].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        }
        if(monomers[nearest].Blist[secBond][1] == third) {
          product(bond, bond0, faceNormal); 
        }
        else {
          product(bond0, bond, faceNormal);
        }
        norm2 = iproduct(faceNormal, faceNormal);
        for(int d=0; d < DIMS; d++)
          faceNormal[d] = faceNormal[d] / sqrt(norm2);

        // 3.2. the relative distance between n1 and the face
        for(int d=0; d < DIMS; d++) {
          v12[d] = monomers[n1].pos_pbc[d] - monomers[nearest].pos_pbc[d];
          if((d==0 && wall_flag < 3) || (d==2 && wall_flag < 2) || (d==1 && wall_flag < 1))
            v12[d] = n_image( v12[d], maxsize[d] );
        }
        dis_beadFace = iproduct(v12, faceNormal);

       // 3.3 L.-J. interaction
        if(dis_beadFace < cutoff)
        {
          if(dis_beadFace < 0.)
            fprintf(stream, "%lf\n", dis_beadFace);
          RepulsiveForce(n1, nearest, monomers, sphere_pm, work_dir);
          RepulsiveForce(n1, second, monomers, sphere_pm, work_dir);
          RepulsiveForce(n1, third, monomers, sphere_pm, work_dir);
        }   
/*
        else {
          for(int d=0; d < DIMS; d++) {
            force[d] = 0.;
            monomers[n1].force[d] += force[d];
            monomers[n1].force_face[d] += force[d];
          }
        }
*/                
      }
      /* 
      else
      {
          /*
          method 1:
          -distance: bead-to-bead
          -direction: the average normal vector of faces around the bead
          -object: n1  and beads in n2's bond list
         
          method 2:
          -object: n1 and beads on the nearest face
          -distance: the distance between n1 and the face
          -direction: the normal vector of the face
          */
    //}
    }            
  }
  fclose(stream);
}

void RepulsiveForce(int n1, int n2, struct monomer *mon, struct sphere_param *sphere_pm, 
     char *work_dir)
{
	extern int max_x, max_y, max_z, wall_flag;
  int maxsize[DIMS];
  double q12[DIMS];
  double q12mag=0.;
	double force[DIMS]={0.};
  double eps= sphere_pm->epsilon_overlap * sphere_pm->kT;;
  double sigma = sphere_pm->eq_overlap / pow(2., 1./6.);
  double cutoff = sphere_pm->eq_overlap;
  double truncation = 0.9*sigma;
  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;
  char filename[100];
  FILE *stream;
  sprintf(filename, "%s/data/dis_beadbead.dat", work_dir);
//  stream = fopen(filename, "a");

	// calculate the monomer-monomer distance
	for(int j=0; j < DIMS; j++) {
	  q12[j] = mon[n1].pos_pbc[j] - mon[n2].pos_pbc[j];
		if((j==0 && wall_flag<3) || (j==2 && wall_flag<2) || (j==1 && wall_flag<1)) 
		  q12[j] = n_image(q12[j], maxsize[j]);     
		q12mag += q12[j]*q12[j];
	}	
	q12mag = sqrt(q12mag);
	// WCA potential
  if(q12mag >= cutoff) { 
    for(int i=0; i<DIMS; i++) 
	    force[i]=0.0; 
  }
	else {
    if(q12mag < truncation)
    {
  	  double aterm = (sigma / truncation);
		  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
		  for(int i=0; i < DIMS; i++) 
		    force[i] = 24.0 * (eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
//fprintf(stream, "truncation=%lf    q12mag=%lf\n", truncation, q12mag);
//printf("force=(%lf,  %lf,  %lf)\n", force[0],force[1],force[2]);
//PAUSE
    }
    else 
    {
	    double aterm = (sigma / q12mag);
		  aterm = aterm*aterm*aterm*aterm*aterm*aterm;
		  for(int i=0; i < DIMS; i++) 
		    force[i]= 24.0 * (eps*(2.0*(aterm*aterm)-aterm)) * q12[i]/q12mag/q12mag;
//fprintf(stream, "q12mag=%lf\n", q12mag);
//printf("force=(%lf,  %lf,  %lf)\n", force[0],force[1],force[2]);
//PAUSE
    }
  }  
	for(int i=0; i < DIMS; i++) {
	  mon[n1].force[i] += force[i];
    mon[n1].force_face[i] += force[i];
    mon[n2].force[i] -= force[i];
    mon[n2].force_face[i] -= force[i];
	}
  for(int i=0; i < DIMS; i++) {
    for(int j=0; j < DIMS; j++) {
      mon[n1].stress_int_v1[i][j] += 0.5*q12[i]*force[j];
      mon[n2].stress_int_v1[i][j] -= 0.5*q12[i]*force[j];
    }
  }        
//  fclose(stream);
}

void ArtificialShift(struct sphere_param *sphere_pm, struct monomer *monomers)
{  
  int totalBead = sphere_pm->num_beads;
  //double sigma = sphere_pm->eq_LJ / pow(2., 1./6.);
  //double cutoff = 1.122 * sigma;
  //double epsWCA = sphere_pm->epsilon_overlap * sphere_pm->kT; 
  //double truncation = 1.0 * sigma;

  for(int n1=0; n1 < totalBead; n1++) 
  { 
    int numNeighbor;
    int nearest;
    double dis_nearest;
    int numBond, secBond, thirdBond, second, third;
    double v12[DIMS];
    double faceNormal[DIMS];
    double faceCOM[DIMS];
    double temp1[DIMS];
    double temp2[DIMS];
    double temp3[DIMS];
    double dis1_inverse=0.;
    double dis2_inverse=0.;
    double dis3_inverse=0.;
    double dis_inverse=0.;
    double dis_beadFace=0.;
    double aterm=0.;  
    double force[DIMS]={0.0};
    extern int max_x, max_y, max_z;
    extern int wall_flag;
    double maxsize[DIMS];
    numNeighbor = monomers[n1].list[0];
    dis_nearest = 100.0;
    maxsize[0]=max_x;
    maxsize[1]=max_y;
    maxsize[2]=max_z;
    int onSameParticle=1;
//printf("numNeighber=%d\n", numNeighbor);

    // 1. find the nearest neighbor, on the other particle, in n1's neighbor list.      
    for(int index=1; index <= numNeighbor; index++) 
    {             
      double dis_temp=0.0;
      int n2 = monomers[n1].list[index];
      if(monomers[n1].sphere_id != monomers[n2].sphere_id) // n1 & n2 are on diff. particles
      {
        onSameParticle=0;  // assign the flag vale

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1)) 
            v12[j] = n_image( v12[j], maxsize[j] );
          dis_temp += v12[j]*v12[j];
        }
        dis_temp = sqrt(dis_temp);

        if(dis_temp < dis_nearest) {
          nearest = n2;
          dis_nearest = dis_temp;
        }
      }
    }// Get the nearest neighbor's index and the distance bwtween them.
       
    // 2. Call nearest neighbor's Bond list to find n1's top 3 nearest neighbors, including n2.    
    if(onSameParticle==0)
    {
      numBond = monomers[nearest].Blist[0][0];
      double dis[numBond+1];
      double dis_temp[numBond];
      int rank[numBond];
      double bond0[DIMS]; 
      double bond[DIMS];
      double norm2=0.;
      dis[0] = dis_nearest;
//printf("numBond=%d\n", numBond);
//printf("dis=(%f %f %f)\n",dis[0],dis[1],dis[2]);
//printf("dis_temp=(%f, %f, %f)\n",dis_temp[0],dis_temp[1],dis_temp[2]);
//printf("rank=(%d %d %d)\n", rank[0],rank[1],rank[2]);

      for(int i=1; i <= numBond; i++)
      {
        int index = monomers[nearest].Blist[i][0];

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[index].pos_pbc[j];  // relative position vector (n1-index)
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
            v12[j] = n_image(v12[j], maxsize[j]);
          dis[i] += v12[j]*v12[j];
        }
        dis[i] = sqrt(dis[i]);
      }
      for(int i=0; i < numBond; i++)
        dis_temp[i] = dis[i+1];

      Sort(dis_temp, numBond, rank);    
      secBond = rank[0];   // the bond index for the second nearest neighbor of n1
      thirdBond = rank[1]; // the bond index for the third nearest neighbor of n1
      second = monomers[nearest].Blist[secBond][0];  // the monomer index of the 2nd nearest neighbor
      third = monomers[nearest].Blist[thirdBond][0];
/*
printf("numBond=%d\n", numBond);
printf("dis[0]=%f dis[1]=%f\n", dis[0],dis[1]);
printf("secBond=%d  thirdBond=%d  second=%d  third=%d\n", secBond, thirdBond, second, third);
PAUSE
*/
      // 3. Calculate the face's normal vector, the relative dist., and the interation
      if(monomers[nearest].Blist[secBond][1] == third || 
         monomers[nearest].Blist[secBond][2] == third) // the three vertexes belong to a face
      {
        // 3.1 normal vector 
        for(int d=0; d < DIMS; d++) {
          bond0[d] = monomers[second].pos_pbc[d] - monomers[nearest].pos_pbc[d];
          bond[d] = monomers[third].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        }
        if(monomers[nearest].Blist[secBond][1] == third) {
          product(bond, bond0, faceNormal); 
        }
        else {
          product(bond0, bond, faceNormal);
        }
        norm2 = iproduct(faceNormal, faceNormal);
        for(int d=0; d < DIMS; d++)
          faceNormal[d] = faceNormal[d] / sqrt(norm2);

        // 3.2. the relative distance between n1 and the face
        for(int d=0; d < DIMS; d++)
          v12[d] = monomers[n1].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        dis_beadFace = iproduct(v12, faceNormal);

        if(dis_beadFace < 0./*0.1*/)
        {
printf("dis_beadFace = %f\n", dis_beadFace);

          //double shift = (0.1-dis_beadFace)/2.;
          //for(int d=0; d < DIMS; d++) {
          //  monomers[n1].pos_pbc[d] += shift * faceNormal[d];
          //  monomers[nearest].pos_pbc[d] -= shift * faceNormal[d]; 
          //  monomers[second].pos_pbc[d] -= shift * faceNormal[d];
          //  monomers[third].pos_pbc[d] -= shift * faceNormal[d];
          //}
          // Modification 20170424 
          for(int d=0; d < DIMS; d++) {       
            monomers[n1].pos_pbc[d] += ((-dis_beadFace)*faceNormal[d]);
            monomers[n1].pos[d] = box(monomers[n1].pos_pbc[d],maxsize[d]);
          } 
        } 
      }
    }            
  }
}

void Artificial_shift(struct sphere_param *sphere_pm, struct monomer *monomers, struct face *faces)
{  
  int totalBead = sphere_pm->num_beads;
  static int counter=0;
  extern char *work_dir;
  int flag=0;

  for(int n1=0; n1 < totalBead; n1++) 
  { 
    int numNeighbor;
    int nearest;
    double dis_nearest;
    int numBond, secBond, thirdBond, second, third;
    double v12[DIMS];
    double faceNormal[DIMS];
    double faceCOM[DIMS];
    double temp1[DIMS];
    double temp2[DIMS];
    double temp3[DIMS];
    double dis1_inverse=0.;
    double dis2_inverse=0.;
    double dis3_inverse=0.;
    double dis_inverse=0.;
    double dis_beadFace=0.;
    double aterm=0.;  
    double force[DIMS]={0.0};
    extern int max_x, max_y, max_z;
    extern int wall_flag;
    double maxsize[DIMS];
    numNeighbor = monomers[n1].list[0];
    dis_nearest = 100.0;
    maxsize[0]=max_x;
    maxsize[1]=max_y;
    maxsize[2]=max_z;
    int onSameParticle=1;

    // 1. find the nearest neighbor, on the other particle, in n1's neighbor list.      
    for(int index=1; index <= numNeighbor; index++) 
    {             
      double dis_temp=0.0;
      int n2 = monomers[n1].list[index];
      if(monomers[n1].sphere_id != monomers[n2].sphere_id) // n1 & n2 are on diff. particles
      {
        onSameParticle=0;  // assign the flag vale

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1)) 
            v12[j] = n_image( v12[j], maxsize[j] );
          dis_temp += v12[j]*v12[j];
        }
        dis_temp = sqrt(dis_temp);

        if(dis_temp < dis_nearest) {
          nearest = n2;
          dis_nearest = dis_temp;
        }
      }
    }// Get the nearest neighbor's index and the distance bwtween them.
       
    // 2. Call nearest neighbor's Bond list to find n1's top 3 nearest neighbors, including n2.    
    if(onSameParticle==0)
    {
      numBond = monomers[nearest].Blist[0][0];
      double dis[numBond+1];
      double dis_temp[numBond];
      int rank[numBond];
      double bond0[DIMS]; 
      double bond[DIMS];
      double norm2=0.;
      dis[0] = dis_nearest;

      for(int i=1; i <= numBond; i++)
      {
        int index = monomers[nearest].Blist[i][0];

        for(int j=0; j < DIMS; j++) {
          v12[j] = monomers[n1].pos_pbc[j] - monomers[index].pos_pbc[j];  // relative position vector (n1-index)
          if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
            v12[j] = n_image(v12[j], maxsize[j]);
          dis[i] += v12[j]*v12[j];
        }
        dis[i] = sqrt(dis[i]);
      }
      for(int i=0; i < numBond; i++)
        dis_temp[i] = dis[i+1];

      Sort(dis_temp, numBond, rank);    
      secBond = rank[0];   // the bond index for the second nearest neighbor of n1
      thirdBond = rank[1]; // the bond index for the third nearest neighbor of n1
      second = monomers[nearest].Blist[secBond][0];  // the monomer index of the 2nd nearest neighbor
      third = monomers[nearest].Blist[thirdBond][0];

      // 3. Calculate the face's normal vector, the relative dist., and the interation

      // ** Note 20170430 ** I can make good use of monomer's data member faceList to 
      // investigate if the three vertices make a triangular face.

      if(monomers[nearest].Blist[secBond][1] == third || 
         monomers[nearest].Blist[secBond][2] == third) // the three vertexes belong to a face
      {
        // 3.1 calculate the normalized normal vector of the face
        for(int d=0; d < DIMS; d++) {
          bond0[d] = monomers[second].pos_pbc[d] - monomers[nearest].pos_pbc[d];
          bond[d] = monomers[third].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        }
        if(monomers[nearest].Blist[secBond][1] == third) 
          product(bond, bond0, faceNormal);         
        else 
          product(bond0, bond, faceNormal);        
        norm2 = iproduct(faceNormal, faceNormal);
        for(int d=0; d < DIMS; d++) 
          faceNormal[d] = faceNormal[d] / sqrt(norm2);  // normalized face normal

        // 3.2. the relative distance between n1 and the face
        // ** Note 20170430 ** Should take the distance of the periodic image into account
        //for(int d=0; d < DIMS; d++)
        //  v12[d] = monomers[n1].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        v12[0] = monomers[n1].pos_pbc[0] - monomers[nearest].pos_pbc[0];          
        v12[1] = monomers[n1].pos_pbc[1] - monomers[nearest].pos_pbc[1];  
        v12[2] = monomers[n1].pos_pbc[2] - monomers[nearest].pos_pbc[2];  
        if(wall_flag < 3)    v12[0] = n_image(v12[0],maxsize[0]);
        if(wall_flag < 2)    v12[2] = n_image(v12[2],maxsize[2]);
        if(wall_flag < 1)    v12[1] = n_image(v12[1],maxsize[1]);
        dis_beadFace = iproduct(v12, faceNormal);

        if(dis_beadFace <= 0.)
        {
          double test[3], a[3], b[3], norm2_a, norm2_b, a_dot_b, a_dot_b_2, test_dot_a, 
                 test_dot_b, denominator, alpha, beta;
          test[0] = v12[0] - dis_beadFace * faceNormal[0];
          test[1] = v12[1] - dis_beadFace * faceNormal[1];
          test[2] = v12[2] - dis_beadFace * faceNormal[2];                          
          //projection[0] = monomers[n1].pos_pbc[0] - dis_beadFace * faceNormal[0];
          //projection[1] = monomers[n1].pos_pbc[1] - dis_beadFace * faceNormal[1];
          //projection[2] = monomers[n1].pos_pbc[2] - dis_beadFace * faceNormal[2];
          //test[0] = projection[0] - monomers[nearest].pos_pbc[0];
          //test[1] = projection[1] - monomers[nearest].pos_pbc[1];
          //test[2] = projection[2] - monomers[nearest].pos_pbc[2];
     
          // ** Note 20170430 ** when determine the basis vectors, do not consider the 
          // periodic image.
          a[0] = monomers[second].pos_pbc[0] - monomers[nearest].pos_pbc[0];
          a[1] = monomers[second].pos_pbc[1] - monomers[nearest].pos_pbc[1];
          a[2] = monomers[second].pos_pbc[2] - monomers[nearest].pos_pbc[2];
          b[0] = monomers[third].pos_pbc[0]  - monomers[nearest].pos_pbc[0];
          b[1] = monomers[third].pos_pbc[1]  - monomers[nearest].pos_pbc[1];
          b[2] = monomers[third].pos_pbc[2]  - monomers[nearest].pos_pbc[2];
          norm2_a = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
          norm2_b = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
          a_dot_b = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
          a_dot_b_2 = a_dot_b * a_dot_b;
          test_dot_a = test[0]*a[0] + test[1]*a[1] + test[2]*a[2];
          test_dot_b = test[0]*b[0] + test[1]*b[1] + test[2]*b[2];
          denominator = norm2_a * norm2_b - a_dot_b_2;      
          alpha = (test_dot_a * norm2_b - test_dot_b * a_dot_b) / denominator;
          beta =  (test_dot_b * norm2_a - test_dot_a * a_dot_b) / denominator;
          if(alpha >= 0. && beta >= 0. && (alpha+beta) <= 1.)
          {
printf("dis_beadFace = %f  alpha=%f  beta=%f\n", dis_beadFace, alpha, beta);
            flag++;
            Write_overlap_config (counter, monomers, faces, sphere_pm, work_dir);
            counter++;
            for(int d=0; d < DIMS; d++) {       
              monomers[n1].pos_pbc[d] += ((-dis_beadFace + 0.1)*faceNormal[d]);
              monomers[n1].pos[d] = box(monomers[n1].pos_pbc[d],maxsize[d]);
            }
            Write_overlap_config (counter, monomers, faces, sphere_pm, work_dir);
            counter++;
          }
          //if(alpha >= 0. && beta >= 0. && (alpha+beta) is much less than 1)
          //{}

        } 
      }
      //else

    }            
  }
}

double rtsafe(/*T &funcd,*/ const double x1, const double x2, const double xacc, 
       struct monomer n1, struct monomer nearest, struct monomer vertex0, 
       struct monomer vertex1) 
{
	const int MAXIT=100;
	double xh,xl;
	double fl = funcd(x1, n1, nearest, vertex0, vertex1);
	double fh = funcd(x2, n1, nearest, vertex0, vertex1);

//extern int n_step;

//	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
//		fatal_err("Root must be bracketed in rtsafe",-1);
////Modification 20170726
//printf("step=%d;  fl=%f;  fh=%f\n",n_step,fl,fh);
//  }
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
//	if (fl < 0.0) {
		xl=x1;
		xh=x2;
//	} else {
//		xh=x1;
//		xl=x2;
//	}
	double rts = 0.5*(x1+x2);
	double dxold = fabs(x2-x1);
	double dx = dxold;
	double f = funcd(rts, n1, nearest, vertex0, vertex1);
	double df = funcd_df(rts, n1, nearest, vertex0, vertex1);
	for (int j=0; j < MAXIT; j++) 
  {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) 
    {
			dxold = dx;
			dx = 0.5*(xh-xl);
			rts =xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			double temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		double f=funcd(rts,n1,nearest,vertex0,vertex1);
		double df=funcd_df(rts,n1,nearest,vertex0,vertex1);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	fatal_err("Maximum number of iterations exceeded in rtsafe",-1);
}

double funcd(double t, struct monomer n1, struct monomer nearest, struct monomer vertex0, 
       struct monomer vertex1)
{
  extern int wall_flag, max_x, max_y, max_z;
  double maxsize[3];
  maxsize[0]=max_x;  maxsize[1]=max_y;  maxsize[2]=max_z;

  double a1[3], d1[3], a2[3], d2[3], r[3], vp1[3];
  for(int n=0; n<3; n++) 
  {
    a1[n]  = vertex1.pos_pbc[n] - nearest.pos_pbc[n];   
    d1[n]  = vertex1.vel[n]     - nearest.vel[n];
    a2[n]  = vertex0.pos_pbc[n] - nearest.pos_pbc[n];
    d2[n]  = vertex0.vel[n]     - nearest.vel[n];
    r[n]   = n1.pos_pbc[n]      - nearest.pos_pbc[n];
    vp1[n] = n1.vel[n]          - nearest.vel[n]; 
  }
  // consider the periodic image
  if(wall_flag < 3)    r[0] = n_image(r[0],maxsize[0]);
  if(wall_flag < 2)    r[2] = n_image(r[2],maxsize[2]);
  if(wall_flag < 1)    r[1] = n_image(r[1],maxsize[1]);  

  double term1[3],term2[3],term3[3],term_a[3],term_b[3];
  product(a1,a2,term1);
  product(a1,d2, term2);
  product(d1,a2, term3);
  for(int n=0; n<3; n++) {
    term2[n] += term3[n];
    term2[n] *= t;
  }
  product(d1,d2,term3);
  term3[0]*=(t*t);  term3[1]*=(t*t);  term3[2]*=(t*t);
  term_a[0]=term1[0]+term2[0]+term3[0];  term_a[1]=term1[1]+term2[1]+term3[1];
  term_a[2]=term1[2]+term2[2]+term3[2];

  vp1[0]*=t;  vp1[1]*=t;  vp1[2]*=t;
  term_b[0]=r[0]+vp1[0];  term_b[1]=r[1]+vp1[1];  term_b[2]=r[2]+vp1[2]; 

  double result = iproduct(term_a,term_b);
  return result;
}

double funcd_df(double t, struct monomer n1, struct monomer nearest, struct monomer vertex0, 
       struct monomer vertex1)
{
  extern int wall_flag, max_x, max_y, max_z;
  double maxsize[3];
  maxsize[0]=max_x;  maxsize[1]=max_y;  maxsize[2]=max_z;

  double a1[3], d1[3], a2[3], d2[3], r[3], vp1[3];
  for(int n=0; n<3; n++) 
  {
    a1[n]  = vertex1.pos_pbc[n] - nearest.pos_pbc[n];   
    d1[n]  = vertex1.vel[n]     - nearest.vel[n];
    a2[n]  = vertex0.pos_pbc[n] - nearest.pos_pbc[n];
    d2[n]  = vertex0.vel[n]     - nearest.vel[n];
    r[n]   = n1.pos_pbc[n]      - nearest.pos_pbc[n];
    vp1[n] = n1.vel[n]          - nearest.vel[n];
  }
  // consider the periodic image
  if(wall_flag < 3)    r[0] = n_image(r[0],maxsize[0]);
  if(wall_flag < 2)    r[2] = n_image(r[2],maxsize[2]);
  if(wall_flag < 1)    r[1] = n_image(r[1],maxsize[1]);  

  double temp1[3],temp2[3];
  product(a1,d2,temp1);
  product(d1,a2,temp2);
  temp1[0]+=temp2[0];  temp1[1]+=temp2[1];  temp1[2]+=temp2[2];
  double term_1 = iproduct(temp1,r);

  product(d1,d2,temp1);
  temp2[0]=r[0]*2*t;  temp2[1]=r[1]*2*t;  temp2[2]=r[2]*2*t;
  double term_2 = iproduct(temp1,temp2);

  product(a1,a2,temp1);
  double term_3 = iproduct(temp1,vp1);

  product(a1,d2,temp1);
  product(d1,a2,temp2);
  temp1[0]+=temp2[0];  temp1[1]+=temp2[1];  temp1[2]+=temp2[2];
  temp2[0]=vp1[0]*2*t;  temp2[1]=vp1[1]*2*t;  temp2[2]=vp1[2]*2*t;
  double term_4 = iproduct(temp1,temp2);

  product(d1,d2,temp1);
  temp2[0]=vp1[0]*3*t*t;  temp2[1]=vp1[1]*3*t*t;  temp2[2]=vp1[2]*3*t*t;
  double term_5 = iproduct(temp1,temp2); 

  double result = term_1 + term_2 + term_3 + term_4 + term_5;
  return result;
}

void preclude_penetraction(struct sphere_param *sphere_pm, struct monomer *monomers, 
     struct face *faces, double dt)
{  
  extern int max_x, max_y, max_z, wall_flag;
  double maxsize[3];
  maxsize[0]=max_x;  maxsize[1]=max_y;  maxsize[2]=max_z;
  char filename[100];
  FILE *stream;
  extern char *work_dir;  

  int totalBead = sphere_pm->num_beads;
  for(int n1=0; n1 < totalBead; n1++) 
  {
    /* Note 20170726: make # of neighbors in a reasonable range
       Might need to modify the code of the neighbor list 
    */
    int breakFlag=FALSE;  // Modification 20170818

    int numNeighbor = monomers[n1].list[0];
    double dis[numNeighbor];
    int label[numNeighbor];
    for(int index=1; index <= numNeighbor; index++) 
    {             
      int n2 = monomers[n1].list[index];
      label[index-1]=n2; 
      double dis_temp=0.0;
      double v12[3];       
      for(int j=0; j < DIMS; j++) {
        v12[j] = monomers[n1].pos_pbc[j] - monomers[n2].pos_pbc[j];
        if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1)) 
          v12[j] = n_image(v12[j], maxsize[j]);
        dis_temp += v12[j]*v12[j];
      }
      dis_temp = sqrt(dis_temp);
      dis[index-1] = dis_temp;      
    }
    sort_label(dis,label);

    int top6=6, temp;       
    if(numNeighbor<top6) temp=numNeighbor;  
    else                 temp=top6;         
    for(int neighbor=0; neighbor<temp; neighbor++)
    {
      if(breakFlag == TRUE) // Modification 20170818
        break;
      else 
      {
        int n2 = label[neighbor];
        if(monomers[n1].sphere_id != monomers[n2].sphere_id)
        {
          // bond label starts from 1
          for(int index=1; index <= monomers[n2].Blist[0][0]; index++)         
          {
            int faceLabel = monomers[n2].faceList[index];
            //n2 = faces[faceLabel].vertices[0]
  //          int vertex0 = faces[faceLabel].vertices[1];
  //          int vertex1 = faces[faceLabel].vertices[2];
  int vertex0 = monomers[n2].Blist[index][0];
  int vertex1 = monomers[n2].Blist[index][1];
            double bond0[3], bond1[3], faceNormal[3], r[3]; 
            bond1[0] = monomers[vertex1].pos_pbc[0] - monomers[n2].pos_pbc[0];
            bond1[1] = monomers[vertex1].pos_pbc[1] - monomers[n2].pos_pbc[1];
            bond1[2] = monomers[vertex1].pos_pbc[2] - monomers[n2].pos_pbc[2];
            bond0[0] = monomers[vertex0].pos_pbc[0] - monomers[n2].pos_pbc[0];
            bond0[1] = monomers[vertex0].pos_pbc[1] - monomers[n2].pos_pbc[1];
            bond0[2] = monomers[vertex0].pos_pbc[2] - monomers[n2].pos_pbc[2];
            product(bond1,bond0,faceNormal);
            r[0] = monomers[n1].pos_pbc[0] - monomers[n2].pos_pbc[0];
            r[1] = monomers[n1].pos_pbc[1] - monomers[n2].pos_pbc[1];
            r[2] = monomers[n1].pos_pbc[2] - monomers[n2].pos_pbc[2];
            
            // Modification 20170818
            //// consider the periodic image
            //if(wall_flag < 3)    r[0] = n_image(r[0],maxsize[0]);
            //if(wall_flag < 2)    r[2] = n_image(r[2],maxsize[2]);
            //if(wall_flag < 1)    r[1] = n_image(r[1],maxsize[1]);
            double b0 = iproduct(r,faceNormal);
  //printf("faceLabel:%d\n",faceLabel);
  //printf("(n2, vertex0, vertex1)=(%d, %d, %d)\n",n2,vertex0,vertex1);
  //printf("n1=(%le, %le, %le)\n",monomers[n1].pos_pbc[0],monomers[n1].pos_pbc[1],monomers[n1].pos_pbc[2]);
  //printf("n2=(%le, %le, %le)\n",monomers[n2].pos_pbc[0],monomers[n2].pos_pbc[1],monomers[n2].pos_pbc[2]);
  //printf("vertex0=(%le, %le, %le)\n",monomers[vertex0].pos_pbc[0],monomers[vertex0].pos_pbc[1],monomers[vertex0].pos_pbc[2]);
  //printf("vertex1=(%le, %le, %le)\n",monomers[vertex1].pos_pbc[0],monomers[vertex1].pos_pbc[1],monomers[vertex1].pos_pbc[2]);
  
            // Vertex positions at the next time step
            double n1_pos_pbc[3], n2_pos_pbc[3], vertex1_pos_pbc[3], vertex0_pos_pbc[3]; 
            for(int d=0; d<3; d++)
            { 
              n1_pos_pbc[d] = monomers[n1].pos_pbc[d] + monomers[n1].vel[d]*dt;
              n2_pos_pbc[d] = monomers[n2].pos_pbc[d] + monomers[n2].vel[d]*dt;
              vertex0_pos_pbc[d] = monomers[vertex0].pos_pbc[d] + monomers[vertex0].vel[d]*dt;
              vertex1_pos_pbc[d] = monomers[vertex1].pos_pbc[d] + monomers[vertex1].vel[d]*dt;
              bond1[d] = vertex1_pos_pbc[d] - n2_pos_pbc[d];
              bond0[d] = vertex0_pos_pbc[d] - n2_pos_pbc[d];
            }
            product(bond1,bond0,faceNormal);
            r[0] = n1_pos_pbc[0] - n2_pos_pbc[0];//monomers[n2].pos_pbc[0];        
            r[1] = n1_pos_pbc[1] - n2_pos_pbc[1];//monomers[n2].pos_pbc[1];        
            r[2] = n1_pos_pbc[2] - n2_pos_pbc[2];//monomers[n2].pos_pbc[2];
            // Modification 20170818 
            //// consider periodic image
            //if(wall_flag < 3)    r[0] = n_image(r[0],maxsize[0]);
            //if(wall_flag < 2)    r[2] = n_image(r[2],maxsize[2]);
            //if(wall_flag < 1)    r[1] = n_image(r[1],maxsize[1]);
            double b1 = iproduct(r,faceNormal);
  
            if(b0*b1<=0.)      
            {
  //printf("b0=%le;  b1=%le\n",b0,b1);
  //PAUSE
  //Modification 20170820 1->1.5
              // predict the time at which the vertex bounds the face 
              double t_prime = rtsafe(0.,1.5,10e-4,monomers[n1],monomers[n2],
                               monomers[vertex0],monomers[vertex1]);
  // Modification 20170726: TEST
  
  //printf("t_prime=%le\n",t_prime);
  
              double g[3];
              for(int d=0; d<3; d++)
              { 
                n1_pos_pbc[d] = monomers[n1].pos_pbc[d] + monomers[n1].vel[d]*t_prime;
                n2_pos_pbc[d] = monomers[n2].pos_pbc[d] + monomers[n2].vel[d]*t_prime;
                vertex0_pos_pbc[d] = monomers[vertex0].pos_pbc[d] + 
                                     monomers[vertex0].vel[d]*t_prime;
                vertex1_pos_pbc[d] = monomers[vertex1].pos_pbc[d] + 
                                     monomers[vertex1].vel[d]*t_prime;
                g[d] = n1_pos_pbc[d]-n2_pos_pbc[d];
              }
              // Modification 20170818 
              // consider the periodic image
              //if(wall_flag < 3)    g[0] = n_image(g[0],maxsize[0]);
              //if(wall_flag < 2)    g[2] = n_image(g[2],maxsize[2]);
              //if(wall_flag < 1)    g[1] = n_image(g[1],maxsize[1]);
  
              // check if n1 is within the face
              double a1[3],a2[3];
              a1[0] = vertex1_pos_pbc[0]-n2_pos_pbc[0];
              a1[1] = vertex1_pos_pbc[1]-n2_pos_pbc[1];
              a1[2] = vertex1_pos_pbc[2]-n2_pos_pbc[2];
              a2[0] = vertex0_pos_pbc[0]-n2_pos_pbc[0];
              a2[1] = vertex0_pos_pbc[1]-n2_pos_pbc[1];
              a2[2] = vertex0_pos_pbc[2]-n2_pos_pbc[2];
              double norm2_a1 = iproduct(a1,a1);
              double norm2_a2 = iproduct(a2,a2);
              double a1_dot_a2 = iproduct(a1,a2);
              double a1_dot_a2_2 = a1_dot_a2 * a1_dot_a2;
              double g_dot_a1 = iproduct(g,a1);
              double g_dot_a2 = iproduct(g,a2);
              double denominator = norm2_a1 * norm2_a2 - a1_dot_a2_2;      
              double zeta = (g_dot_a1 * norm2_a2 - g_dot_a2 * a1_dot_a2) / denominator;
              double xi =   (g_dot_a2 * norm2_a1 - g_dot_a1 * a1_dot_a2) / denominator; 
              if(zeta >= 0. && xi >= 0. && (zeta+xi) <= 1.)
              {
                // update n1 position
                double v_prime[3];  
                for(int d=0; d<3; d++)
                {
                  v_prime[d] = (1-zeta-xi)*monomers[n2].vel[d] + 
                               zeta*monomers[vertex1].vel[d] + xi*monomers[vertex0].vel[d];
                  monomers[n1].vel_temp[d] = 2*v_prime[d] - monomers[n1].vel[d];
                  monomers[n1].pos_temp[d] = n1_pos_pbc[d] + (dt-t_prime) *
                                             monomers[n1].vel_temp[d];
                } 
                // label n1
                monomers[n1].updatedFlag=TRUE;
  
  //printf("t_prime=%le\n",t_prime);
  //PAUSE 
                // Modification 20170818
                breakFlag=TRUE;
                sprintf(filename,"%s/data/t_prime.dat",work_dir);
                stream=fopen(filename,"a");
                fprintf(stream, "%f\n", t_prime);
                fclose(stream);
                break;
              }
            }  
          }        
        }
      }
    }
  }
}

