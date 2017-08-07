#include "header.h"

void CalcuScatIntensity (int size_qx, int size_qy, int size_qz, struct vector *q_vector, 
                         int **node_map, double *intensity, double *i_xz, int output_index, 
                         char *work_dir)
{
  struct vector *marker = NULL;
  int size_node = 0;
  
  marker = GetScatteringSource(node_map, &size_node);
  
  ScatteringIntensity(size_qx, size_qy, size_qz, size_node, q_vector, marker, intensity);

  //SuperimposeIntensity(size_qx, size_qy, size_qz, intensity, i_xz);

  // OutputData(output_index, size_qx, size_qy, size_qz, intensity, i_xz, work_dir);
  
  // free(marker);
}

struct vector * GetScatteringSource (int **node_map, int *count)
{
  extern int max_x, max_y, max_z,  num_x;
  int size_z = max_z+2;
  int size_y = max_y+2;
  int size_x = num_x+2;
  
  struct vector *temp_ptr;
  struct vector *r = NULL;
  *count=0;
  
  for(int k=0; k < size_z; k++)
  {
    for(int j=0; j < size_y; j++)
    {
      for(int i=0; i < size_x; i++)
      {
	int ij = i * size_y + j;

	if (node_map[ij][k]==2)          
	{	
	  (*count)++;
          
          temp_ptr = (struct vector*) realloc (r, (*count)*sizeof(struct vector));
	  
	  if(temp_ptr != NULL)
	  {
	    r = temp_ptr;
	    r[*count-1].x = i;
	    r[*count-1].y = j;
	    r[*count-1].z = k;
	  }
	  
	  else
	  {
	    free(r);
	    printf("Memory allocation for scattering sources failed\n");
	    exit(1);
	  }
	}
      } 		  
    }
  }

  return r;
}

void ScatteringIntensity (int size_qx, int size_qy, int size_qz, int size_node, struct vector *q, 
                          struct vector *r, double *intensity)
{  
  int nxy = size_qx*size_qy;
 
  for(int k=0; k<size_qz; k++)
  {
    for(int j=0; j<size_qy; j++)
    {
      for(int i=0; i<size_qx; i++)
      {
	int index = k*nxy + j*size_qx + i;
	double real_part = 0.0;
	double imaginary_part = 0.0;
        
	for(int n=0; n < size_node; n++)
	{
	  double qr = (q[index].x)*(r[n].x) + (q[index].y)*(r[n].y) + (q[index].z)*(r[n].z); 
	  
	  //printf("r = (%f, %f, %f)\n", r[n].x, r[n].y, r[n].z);
	  real_part +=  sin(qr);
	  imaginary_part += sin(qr);
	}
	
	intensity[index] = (real_part*real_part + imaginary_part*imaginary_part);
	//printf("intensity[index]= %lf\n", intensity[index]);
      }
    }
  }  
}

void SuperimposeIntensity (int size_qx, int size_qy, int size_qz, 
			      double *intensity, double *i_xz)
{
  int nxy = size_qx * size_qy;

  for(int k=0; k<size_qz; k++)
    for(int i=0; i<size_qx; i++)
      for(int j=0; j<size_qy; j++)
      {
	int n = k*size_qx + i;
	int index = k*nxy + j*size_qx + i;
	i_xz[n] += intensity[index];
      }
/*  
  for(int k=0; k<size_qz; k++)
    for(int j=0; j<size_qy; j++)
      for(int i=0; i<size_qx;i++)
      {
	int n = k*size_qy + j;
	int index = k*nxy + j*size_qx + i;
	i_yz[n] += intensity[index];
      }

  for(int j=0; j<size_qy; j++)
    for(int i=0; i<size_qx; i++)
      for(int k=0; k<size_qz; k++)
      {
	int n = j*size_qx + i;
	int index =  k*nxy + j*size_qx + i;

	i_xy[n] += intensity[index];
      }
*/
}

void OutputData (int output_index, int size_qx, int size_qy, int size_qz, 
		 double *intensity, double *i_xz, char *work_dir)
{
  int nxy = size_qx*size_qy;
  FILE *stream;

  for(int j=0; j<size_qy; j++)
  {
    char filename1[500];
    sprintf(filename1, "%s/data/i_y%d_t%d.dat", 
	                          work_dir, j, output_index);
    stream = fopen(filename1, "w");

    for(int k=0; k<size_qz; k++)
    {
      for(int i=0; i<size_qx; i++)
      {
	int index = k*nxy + j*size_qx + i;

        fprintf(stream, "%lf ", intensity[index]);
      }

      fprintf(stream,"\n");
    }

    fclose(stream);
  }

  char filename2[500];
  sprintf(filename2, "%s/data/ixz_%d.dat", 
                                work_dir, output_index);
  stream = fopen(filename2, "w");

  for(int k=0; k<size_qz; k++)
  {
    for(int i=0; i<size_qx; i++)
    {
      int index = k*size_qx + i;
      fprintf(stream, "%lf ", i_xz[index]);
    }
    
    fprintf(stream, "\n");
  }
  
  fclose(stream);  

}

//void GetScatteringSource (int output_index, char *work_dir, struct vector *r)
/*{
  char filename1[500];
  sprintf(filename1, "%s/data/total_nodes_inside_%d.dat", 
	                         work_dir, output_index);
  
  FILE *stream;
  stream = fopen(filename1, "r");
  
  int size;
  fscanf(stream, "%d", &size);
  fclose(stream);
  
  r = (struct vector*)malloc(size*sizeof(struct vector));
  if(r ==NULL)
  {
    printf("Memory allocation for scattering sources failed\n");
    exit(1);
  }
      
  char filename2[500];
  sprintf(filename2, "%s/data/mark_particle_%d.dat", work_dir, output_index);
*/
/*
  stream = fopen(filename2,"r");
  
  for(int n=0; n<size; n++)
  {
    fscanf(stream, "%lf%lf%lf", &r[n].x, &r[n].y, &r[n].z);
  }
  
  fclose(stream);
  free(r);
}
*/
