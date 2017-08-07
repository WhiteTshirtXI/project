# include "header.h"

void AggreForce(struct sphere_param *sphere_pm, struct monomer *monomers)
{  
  int totalBead = sphere_pm->num_beads;
  double sigma = sphere_pm->range_LJ;
  double cutoff = 2.5*sigma;
  double epsLJ = sphere_pm->eps_LJ * sphere_pm->kT; 
 
  for(int n1=0; n1 < totalBead; n1++) 
  { 
    int numNeighbor;
    int nearest;
    double dis_nearest;
    int numBond, secBond, thirdBond, second, third;
    double v12[DIMS];
    double faceNormal[DIMS];
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

    // 1. find the nearest neighbor, on the other particle, in n1's neighbor list.      
    for(int index=1; index <= numNeighbor; index++) 
    {             
      double dis_temp=0.0;
      int n2 = monomers[n1].list[index];
      if(monomers[n1].sphere_id != monomers[n2].sphere_id) // on different particles
      {
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
    numBond = monomers[nearest].Blist[0][0];
    double dis[numBond+1] = {0.0};
    double dis_temp[numBond] = {0.0};
    int rank[numBond];
    double bond0[DIMS]; 
    double bond[DIMS];
    double norm2=0.;
    dis[0] = dis_nearest;

    for(int i=1; i <= numBond; i++)
    {
      int index = monomers[nearest].Blist[i][0];

      for(int j=0; j < DIMS; j++) {
        v12[j] = monomers[n1].pos_pbc[j] - monomers[index].pos_pbc[j];  // n1-index
        if((j==0 && wall_flag < 3) || (j==2 && wall_flag < 2) || (j==1 && wall_flag < 1))
          v12[j] = n_image(v12[j], maxsize[j]);
        dis[i] += v12[j]*v12[j];
      }
      dis[i] = sqrt(dis[i]);
    }
    for(int i=0; i < numBond; bond++)
      dis_temp[i] = dis[i+1];

    Sort(dis_temp, numBond, rank);    
    secBond = rank[0];   // the bond index for the second nearest neighbor
    thirdBond = rank[1]; // the bond index for the third nearest neighbor
    second = monomers[nearest].Blist[secondBond][0];  // monomer indexes
    third = monomers[nearest].Blist[thirdBond][0];

    // 3. Calculate the face's normal vector, the relative dist., and the interation
    if(monomers[nearest].Blist[secondBond][1] == third || 
       monomers[nearest].Blist[secondBond]][2] == third) // the three vertexes belong to a face
    {
      // 3.1 normal vector 
      for(int d=0; d < DIMS; d++) {
        bond0[d] = monomers[second].pos_pbc[d] - monomers[nearest].pos_pbc[d];
        bond[d] = monomers[third].pos_pbc[d] - monomers[nearest].pos_pbc[d];
      }
      if(monomers[nearest].Blist[secondBond][1] == third) {
        product(bond, bond0, faceNormal); 
      }
      else {
        product(bond0, bond, faceNormal);
      }
      norm2 = iproduct(faceNormal, faceNormal);
      faceNormal[d] = faceNormal[d] / sqrt(norm2);

      // 3.2. the relative distance between n1 and the face
      for(int d=0; d < DIMS; d++)
        v12[d] = monomers[n1].pos_pbc[d] - monomers[nearest].pos_pbc[d];
      dis_beadFace = iproduct(v12, faceNormal);

      // 3.3 L.-J. interaction
      if(dis_beadFace < cutoff)
      {
        // Todo: n1 receives a force and the vorteices of the face experience a reaction force.        
         aterm = (sigma / dis_beadFace);
         aterm = aterm * aterm * aterm * aterm * aterm * aterm;
         for(int d = 0; d < DIMS; d++)
           force[d] = 24.0 * (epsLJ/dis_beadFace) * (2.0*(aterm*aterm) - aterm) * faceNormal[d];
         for(int d=0; d < DIMS; d++)
           monomers[n1].force[d] += force[d];
      }          
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
