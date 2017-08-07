 /*  Structure definition header file  */
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


struct vector  {double  x, y, z;};
struct smatrix {double  xx, yy, zz, yz, zx, xy;};
struct matrix  {double  xx, xy, xz, yx, yy, yz, zx, zy, zz;};

struct object                                    /* Objects */
{
  int (*func_ptr)
      (struct object *, double, double, double);   /* Defining function */

  struct vector  r, e;                             /* Coordinates */
  struct vector  r_map, r_lst;                     /* Extra coordinates */
  struct vector  rold;                             /* t-1 coordinates */
  struct vector  u, w;                             /* Velocities */
  struct vector  f, t;                             /* Forces */
  struct vector  p, l;                             /* Accumulated momentum */
  struct vector  f_ext;                            /* External forces */
  struct vector  c, rxc;                           /* Moments of boundary nodes */

  struct smatrix  ztt, zrr, mtt, mrr;              /* Moments and inverse */
  struct matrix   ztr, mtr;                        /* Moments and inverse */
  struct smatrix pf_str, pf_lub;                   /* Particle-fluid stress */
  struct smatrix pc_str, p_str;                    /* Particle-particle stress */
  struct e_diss {double t, r, c;} e_diss;          /* Energy dissipation */

  double mass, inertia;                            /* Mass, Inertia */
  double r_a, r_b;                                 /* Major and Minor Axes */
  double rho;                                      /* Particle density ratio */
  double lub_cut;                                  /* Lubrication cut off */
  double del_hy;                                   /* Hydrodynamic layer */
  double vol;                                      /* Internal volume */
  double sum;                                      /* Sum of fluid bnodes */
  double fluid;                                    /* Fluid mass inside particle */
  double dr2;                                      /* mean square displacement (r(t)-r(t-1))^2 */
  double disp2;                                    /* mean squared displacement (r-r(0))^2 */

  int    shape_flag;                               /* 0=sph 1=ellip 2=open cyl 3=closed cyl 4=capped cyl */
  int    mass_flag;                                /* 0 = Infinite mass */
  int    move_flag;                                /* 0 = Stationary objects */
  int    max_bnode;                                /* Max # of nodes */
  int    size;                                     /* Max object dimension */
  int    list[MAX_N+1];                            /* Neighbor list */
  int    cl_list[MAX_N+1];                         /* Linked particle list */
  int    n_cl;                                     /* Cluster index */
  int    *nodes, *links, num_bnode;                /* Node coords, links, # boundary nodes */
  int    num_coll;                                 /* # of collisions in last Num_Lbe_Step */
  int    stop_flag;                                /* Stopped particle flag (1 = stopped) */
  struct vector u_init;                            // ctliao 20151110
};

struct list                                      /* Collision list */
{
  double t_coll;
  int    coll_1, coll_2, next;
};

struct cluster                                   /* Cluster list */
{
  int *list;                                       /* Cluster list pointers */
};

struct monomer                                   /* monomers */
{
  long sphere_id;                                /* sphere id*/
  int type;                                      /* sphere type */
  double radius;                                 /* monomer radius */
  double rho;                                    /* monomer density */
  double pos[DIMS];                              /* monomer position */
  double pos0[DIMS];                             /* monomer position */
  double pos_pbc[DIMS];                          /* pbc monomer pos */
  double pos_lst[DIMS];                          /* pbc monomer pos */
  double pos_old[DIMS];                          /* pbc monomer pos */
  double pos_tmp[DIMS];                          /* pbc monomer pos */
  double dr2;                                    /* pbc monomer displacement */

  double vel[DIMS];                              /* monomer velocity */
  double force[DIMS];                            /* monomer force */
  double force0[DIMS];
  double fricforce[DIMS];                        /* friction force on the monomer */
  double fluid_vel[DIMS];                        /* fluid velocity at the monomer position */
  int    list[MAX_N+1];                          /* Neighbor list */
  int    blist[MAX_BOND+1][3];              /* bond list (first element is the bonded monomer; second and third are faces) */
  double initLength[MAX_BOND+1];  // initial bond length. 1st element contains nothing. 20170219
  double initAngle[MAX_BOND+1];   // initial bending angle. 20170219 
  double initLength_temp[MAX_BOND+1]; // 20170308
  double initAngle_temp[MAX_BOND+1];  // 20170308
  double initLength_final[MAX_BOND+1];
  double f_ext[DIMS];                            /* extern force acting on the monomers */
  double f_int[DIMS];                            /* internal forces */
  double f_fluc[DIMS];                           /* fluctuation forces */
  double stress[DIMS][DIMS];
  double force_pp[3];
  double drag[DIMS];
  
  int Blist[MAX_BOND+1][3]; // all the bonds of a vertex
  int faceList[MAX_BOND+1];
  int num_proj;
  double Dis[MAX_F+1];
  double f_aggre[DIMS];
  double curve;
  double normal[DIMS];
  double force_spring[DIMS];
  double force_vol[DIMS];
  double force_bend[DIMS];
  double force_area[DIMS];
  double force_ev[DIMS];
  double force_wall[DIMS];
  double force_fluc[DIMS];
  double force_face[DIMS];
  double stress_elas[3][3];
  double stress_bend[3][3];
  double stress_vol[3][3];
  double stress_area[3][3];
  double stress_wall[3][3];
  double stress_int_v1[3][3];
  double stress_int_v2[3][3];
  double dr; // 20161129
  double v_diff[3]; // 20170125
  double pos_com[DIMS]; // 20170220
  double force_nonbonded_intra[3];
  double vel_old[3];
  int updatedFlag;  // Modification 20170725
};

struct face
{
  int sphere_id;                                   /* face id*/
  int vertices[3];                            /* vertices of the face */
  double area;                                   /* area of the face */
  double area_0;
};

struct sphere                                    /* spheres */
{
  double com[DIMS], com0[DIMS], comold[DIMS];    /* center of mass */
  double Rg2;                                    /* radius of gyration */
  double Rx2, Ry2, Rz2;                          /* principal moments of gyration tensor */
  double Ia, Ib, Ic;                             /* principal moments of inertia */
  double asphericity;
  double acylindricity;
  double stretch[DIMS];                          /* longest distance between two beads */
  double avg_sprng_len;                          /* average spring length */
  double disp2, dr2;
  double theta;                                  /* angle between the eigenvector of inertia tensor and flow direction */
  double omega[DIMS];                            /* angular velocity */
  double area;                                   /* total surface area (ie, sum of triangle areas) */
  double volume;                                 /* total volume, calculated by sum of volumes of triangular cones */
  double volume_dupin;                           /* total volume, calculated using Dupin's method */
  double tcom[DIMS], tcomold[DIMS];              /* temp center of mass for node update */
  double tstretch[DIMS], tstretch_old[DIMS];     /* temp stretch for node update */
  
  double g_eigenvalue[DIMS];
  double i_eigenvector[9];
  double semiaxis_a;
  double semiaxis_b;
  double semiaxis_c;
  double deformation; // 20170131
  double inclination; // 20170131
};

struct sphere_param                              /* parameters for sphere */
{
  int spring;                                    /* spring type */
  int ev_type;                                   /* exc. vol. type */
  int verlet;                                    /* verlet update type */
  int Nsphere;                                   /* # of spheres */
  int Ntype[NTYPES];                             /* # of spheres for each type */
  int nlevel[NTYPES];                            /* level of refinement in triangulation (0=icosahedra) for each type */
  int N_per_sphere[NTYPES];                      /* # of beads per sphere */
  int face_per_sphere[NTYPES];                   /* # of faces per sphere */
  int num_beads;                                 /* total # of beads */
  int nfaces;                                    /* total # of faces */
  int numsteps;                                  /* total # of LBE steps */
  int MD_steps;                                  /* # of MD steps per LBE step */
  int initconfig;                                /* sphere initial configuration */
  int relax_time;                                /* sphere relaxation time */
  int write_time;                                /* timestep between output properties */
  int write_config;                              /* timestep between output configs */
  int write_fluid;                               /* timestep between output fluid */

  
  double sigma_k;                                /* Kuhn segment length in units of lattice spacing */
  double H_fene[NTYPES], Q_fene[NTYPES];         /* FENE parameters */
  double k_bend[NTYPES];                         /* bending force strength */
  double k_V[NTYPES];                            /* volume restoring force strength */
  double k_A[NTYPES];                            /* area restoring force strength */
  double Ss;                                     /* Ideal chain bead radius */

  double monmass;                                /* monomer mass */
  double fric;                                   /* friction coef. on a point force */
  double evcutoff;                               /* cutoff for exc. vol. interactions */
  double nks;                                    /* # of kuhn seg. / spring */
  double dt;                                     /* verlet timestep */
  double kT;                                     /* temperature */
  double tempscale;                              /* ratio between input and output monomer temperature */
  double f_ext[DIMS];                            /* extern force acting on the monomers */
  double V0[NTYPES];                             /* rest state volume */
  double A0[NTYPES];                             /* rest state area */
  double V0_temp[NTYPES]; // 20170308
  double A0_temp[NTYPES]; // 20170308 
  double V0_final[NTYPES];
  double A0_final[NTYPES];
  double V0_init[NTYPES];                        // ctliao 20151110
  double cutoff_LJ;
  double epsilon_LJ;
  double cutoff_LJ_final;
  double epsilon_LJ_final;
//  double range_LJ;
  double eq_LJ;
  double epsilon_WCA;
  double eq_WCA;          // Modification 20170725
  double epsilon_overlap;
  double eq_overlap;
  int attrac_type;
  double eq_morse;
  double depth_morse;
  double width_morse;
  double cutoff_morse;  
  double persisLength; 
  double ka_local;
  int attrac_type_1;  // Modification 20170725
};

