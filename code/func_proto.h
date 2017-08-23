/*  Function prototype header file  */
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
#include "mkl_vsl.h"

void   driver (char *);

int    update (int, int, int, int, int, int, int, int, VSLStreamStatePtr, int **, Float ***, 
               double [MAX_Y][Num_Prop], char *, char *, char *, char *,  struct vector, 
               struct object *, struct sphere *, struct monomer *, struct face *, 
               struct sphere_param *, char *, int, int writeInterval_stress, int retrieve_flow_step, struct vector **forceDen);  // Modification 20170419

void lbe_update(struct object *, struct sphere_param *, struct sphere *, struct monomer *, 
                Float ***, int **, struct vector, struct vector *, double [MAX_Y][16], int, 
                int, FILE *, VSLStreamStatePtr, double [6], double [6]);
Float aWeight(int idVel);
void   spheroid_init (struct object *, double, double);
void   wall_init     (struct object *, int);
int    sphere     (struct object *, double, double, double);
int    ellipsoid  (struct object *, double, double, double);
int    cylinder   (struct object *, double, double, double);
int    capped_cyl (struct object *, double, double, double);
int    wall       (struct object *, double, double, double);
int    sphere_init(struct sphere_param *, struct sphere *, struct monomer *, struct face *, char *);
int    check_overlap(struct monomer, struct monomer);
int    check_walloverlap(struct monomer) ;
void   get_forces(struct sphere_param *, struct monomer *, struct face *, Float ***, int, VSLStreamStatePtr);
void   verlet_update(struct monomer *, struct face *, struct sphere_param *, Float ***, int, VSLStreamStatePtr);
void   vel_fluc(struct monomer *);
void   sphere_props(struct sphere_param *, struct sphere *, struct monomer *, struct face *);
void   sphere_props_short(struct sphere_param *, struct sphere *, struct monomer *);
void   bnodes_init  (struct object *, int **);
void   bnodes_add   (struct object *, Float ***);
void   bnodes_del   (struct object *, Float ***);
void   bnodes_mom   (struct object *, int **);
void   bnodes_sph_1 (struct object *, Float ***, int **);
void   bnodes_sph_2 (struct object *, Float ***, int **);
void   bnodes_wall  (struct object *, Float ***, int **);
void   bnodes_dp    (int **, struct sphere_param *, struct sphere *, struct monomer *);
void   implicit_force (struct object *, struct vector);
void   matrix_inv     (struct object *);
void   lbe_bconds  (Float ***);
void   lbe_move  (Float ***);
void lbe_zcol(Float **, int *, struct vector, double [16], int, double [], double [6], 
              double [6]);
void   modes_write (Float **, int *, struct vector, FILE *);
void   check_vel_conservation(Float ***, int);
void   n_list (struct object *);
void   n_list_mon (struct sphere_param *sphere_pm, struct monomer *mon);
void   lub    (struct object *, double);
void   velcs_update  (struct object *, double);
void   hs3d   (struct object *, double);
int    coll   (struct object *, struct list *, double *, double, double, int, int, int, int);
void   cluster_index  (struct object *, int *);
void   cluster_make   (struct object *, struct cluster *, int);
void   cluster_force  (struct object *, struct cluster *, double, int, int);
void   cluster_update (struct object *, struct cluster *, double, int);
void   cluster_matrix (struct object *, struct cluster *, double **, int, double);
void   cj_grad   (double **, double *, int);
void   globals   (struct object *);
void   force_sum (struct object *);
void   multi_sum (double *, int);
void   file_name (char *, char *, int);
void   gauss_jordan (double [6][12]);
void   warning   (char *);
void   fatal_err (char *, int);
void   error_chk ();

void  Write_Output(int, int, int, int, int, int, int, struct object *, struct sphere *, 
                   struct monomer *, struct face *, struct sphere_param *, Float ***, int **, char *);

void   Write_Particle(struct object *, struct sphere_param *, int, char *);
void   Write_Fluid(Float ***, int **, int, char *);
void   Write_Monomer(struct sphere *, struct monomer *, struct sphere_param *, int , int, char * );
void   Write_Sphere(struct sphere *, struct monomer *, struct face *, struct sphere_param *, int , int, char * );
void   Write_Sphereconfig(struct sphere *, struct monomer *, struct face *, struct sphere_param *, int , int , char *);
void   Write_stat(Float ***, int, double, char *);
void   Write_Velslice(struct monomer *, Float ***, int **, int, char *);
void   write_velocity_field(int, Float ***, int **, char *);
void   Write_Nodemap(int **, int, char *);
void   Write_time(int, char *);
void   Write_RDF(struct sphere * , struct monomer *, struct sphere_param * , int , int, int, char *);
void   vector_copy (double *, int, int, int, int);
void   vector_xchg (double *, double *, int, int, int);
void   broad_cast  (double *, int, int);
void   global_sum  (double *, int);
int    global_max  (int);
void   init_procs  (int*, char ***);
void   fini_procs  ();
void   sync_procs  ();
int    proc_num ();
int    proc_id  ();
double wclock ();
double rand_num(long *, double,double);
void SetFluctuationAmplitudes(Float *, Float, Float, Float, Float);
void GenerateFluctuations(Float *);
void CheckVslError(int);
void product(double a[DIMS], double b[DIMS], double c[DIMS]);

void  MarkParticle (int, int **, struct sphere_param*, struct monomer*, char *);

void  CalcuScatIntensity (int, int, int, struct vector *, int **, double *, double *, int, char *);

struct vector *GetScatteringSource (int **, int *);

void  ScatteringIntensity (int, int, int, int, struct vector *, struct vector *, double *);

void  SuperimposeIntensity (int, int, int, double *, double *);

void  OutputData (int, int, int, int, double *, double *, char *);

void Discretize (int flow_flag, int window, int mark_interval, int oscillation_period, 
		 int time_step, int **node_map, struct sphere_param *sphere_pm, 
		 struct monomer *monomers, char *work_dir);
void WallStress(int step, struct vector wallVel, Float ***velcs_df, char *work_dir, int
                outputInterval);

void WriteParticleStress(int step, int outputInterval, struct sphere_param *sphere_pm,
     struct monomer *monomers, char *work_dir);

void WriteBlist(struct sphere_param *sphere_pm, struct monomer *monomers, char *work_dir);
void AssignBlist(struct sphere_param *sphere_pm, struct monomer *monomers, char *work_dir);
void SetFace(struct sphere_param *sphere_pm, struct monomer *monomers, struct face *faces);
double iproduct(double a[DIMS], double b[DIMS]);
void AggreForce(struct sphere_param *sphere_pm, struct monomer *monomers);
int* Sort(double arr[], int len, int *index);
void WriteForce(int step, struct sphere_param *sphere_pm, struct monomer *mono, char *work_dir);
void AggreForce2(struct sphere_param *sphere_pm, struct monomer *monomers, char *work_dir);
void RepulsiveForce(int n1, int n2, struct monomer *mon, struct sphere_param *sphere_pm,
     char *work_dir);
void ArtificialShift(struct sphere_param *sphere_pm, struct monomer *monomers);
void WriteTemplate(struct sphere_param *sphere_pm, struct monomer *mono, char *work_dir);
void ParticleDim(struct sphere_param *sphere_pm, struct monomer *mono);
void BubbleSort(double *arr, int len);
void write_particle_para (struct sphere_param *sphere_pm, struct monomer *vertex, struct face *faces, char *work_dir);
void biconcaveTemplate (struct sphere_param *sphere_pm, struct monomer *vertex);
void read_particle_para (struct sphere_param *sphere_pm, struct monomer *vertex, struct face *faces, char *work_dir);
void calculate_particle_para (struct sphere_param *sphere_pm, struct monomer *vertex, struct face *faces, char *work_dir);
void nonbonded_interaction(struct sphere_param *sphere_pm, struct monomer *monomers);
void nonbonded_interaction_nlist(struct sphere_param *sphere_pm, struct monomer *monomers);
//double* distribution_eq(int xy, int z, double ***velcs_df);
void spreading           (struct vector **forceDen, struct monomer *mon, int num_beads);
void interpolation       (int num_beads, Float ***velcs_df, struct vector **forceDen, double dt, struct monomer *mon);
void equilibrium_distrib (int xy, int z, double ***velcs_df, double dt, struct vector forceDen, struct vector *correctedVel, double *f_eq);
void external_force      (double tau, struct vector forceDen, struct vector correctedVel, double *f_ext);
void collision (double tau, double ***velcs_df, struct vector **forceDen, double dt, int step, int writeInterval);
double adams_bashforth   (double y, double derivOld, double derivNew, double dt);
double euler_method      (double y, double deriv, double dt);
void update_position     (int numBead, double dt, struct monomer *mon);
void propagation         (struct object *objects, int ** node_map, Float ***velcs_df);
void growth_procedure    (int n_step, int totalGrowthStep, struct sphere_param *sphere_pm,
                          struct monomer *monomers, struct object *objects, int n_cycle, 
                          int num_step, int mark_interval, int oscillation_period,  
                          int window, struct sphere *spheres, struct face *faces, 
                          Float ***velcs_df, int **node_map, char *work_dir, VSLStreamStatePtr rngstream);
void euler_update (int numBead, struct monomer *mon);
void growth_procedure_2 (int n_step, int totalGrowthStep, struct sphere_param *sphere_pm, struct monomer *monomers, 
                         struct object *objects, int n_cycle, int num_step, int mark_interval, int oscillation_period, 
                         int window, struct sphere *spheres, struct face *faces, Float ***velcs_df, int **node_map, 
                         char *work_dir, VSLStreamStatePtr rngstream);
void get_forces_hi (struct sphere_param *, struct monomer *, struct face *, Float ***, int, VSLStreamStatePtr);
void Write_overlap_config(int n,struct monomer *monomers, struct face *faces, struct sphere_param *sphere_pm, char *work_dir);
void Artificial_shift    (struct sphere_param *sphere_pm, struct monomer *monomers, struct face *faces);
void write_growth_config (struct sphere *spheres,struct monomer *monomers,struct face *faces,struct sphere_param *sphere_pm,int n_step, 
                          int addsteps,	char *work_dir);
void get_force_growth (struct sphere_param *sphere_pm,struct monomer *monomers,struct face *faces,Float ***velcs_df,int n_step, 
                       VSLStreamStatePtr rngstream);
double funcd(double t, struct monomer n1, struct monomer nearest, struct monomer vertex0, struct monomer vertex1);
double funcd_df(double t, struct monomer n1, struct monomer nearest, struct monomer vertex0, struct monomer vertex1);
double rtsafe(/*T &funcd,*/ const double x1, const double x2, const double xacc, struct monomer n1, struct monomer nearest, 
               struct monomer vertex0, struct monomer vertex1);
void preclude_penetraction(struct sphere_param *sphere_pm, struct monomer *monomers, struct face *faces, double dt);
void sort_label(double arr[], int *index);

