#include "header.h"

void face_props(struct sphere_param *sphere_pm, struct monomer *mons, struct face *faces)
{
	extern int max_x, max_y, max_z;
	extern int wall_flag;
	extern double maxsize[DIMS];
	int i, j, k, d, n;
	int n1, n2, n3, f2;
	int num_beads = sphere_pm->num_beads;
	int Nsphere[NTYPES];
	int Nbead, Nface, nface, numsphere, type, face0, bead0;
	double dr[DIMS];
	double normal[DIMS], normal_mag;
	double **com, *delta_V, *delta_A;
	double q1[DIMS], q2[DIMS];

	for (j = 0; j < NTYPES; j++)
		Nsphere[j] = sphere_pm->Ntype[j];
	numsphere = Nsphere[0] + Nsphere[1];

	com = (double **) calloc(numsphere, sizeof(double *));
	delta_V = (double *) calloc(numsphere, sizeof(double));
	delta_A = (double *) calloc(numsphere, sizeof(double));
	for (i = 0; i < numsphere; i++)
		com[i] = (double *) calloc(DIMS, sizeof(double));

	for (i = 0; i < numsphere; i++) {
		type = (i < sphere_pm->Ntype[0] ? 0 : 1);
		Nbead = sphere_pm->N_per_sphere[type];
		Nface = sphere_pm->face_per_sphere[type];

		face0 = (type == 0 ? i * sphere_pm->face_per_sphere[0] : Nsphere[0] * sphere_pm->face_per_sphere[0] + (i - Nsphere[0]) * sphere_pm->face_per_sphere[1]);
		bead0 = (type == 0 ? i * sphere_pm->N_per_sphere[0] : Nsphere[0] * sphere_pm->N_per_sphere[0] + (i - Nsphere[0]) * sphere_pm->N_per_sphere[1]);

		/* calculate the center-of-mass of each sphere */
		for (d = 0; d < DIMS; d++) {
			com[i][d] = 0.0;
			for (j = 0; j < Nbead; j++)
				com[i][d] += mons[bead0 + j].pos_pbc[d];
			com[i][d] /= Nbead;
			for (j = 0; j < Nbead; j++)	// 20141028 ctliao
				mons[bead0 + j].posCM[d] = mons[bead0 + j].pos_pbc[d] - com[i][d];
		}

		/* calculate face area, total volume and normal vector */
		delta_V[i] = 0.;
		delta_A[i] = 0.;
		for (j = 0; j < sphere_pm->face_per_sphere[type]; j++) {
			nface = face0 + j;
			n1 = faces[nface].vertices[0];
			n2 = faces[nface].vertices[1];
			n3 = faces[nface].vertices[2];

			for (d = 0; d < DIMS; d++)
				dr[d] = mons[n1].pos_pbc[d] - com[i][d];
			for (d = 0; d < DIMS; d++)
				q1[d] = mons[n2].pos_pbc[d] - mons[n1].pos_pbc[d];
			for (d = 0; d < DIMS; d++)
				q2[d] = mons[n3].pos_pbc[d] - mons[n1].pos_pbc[d];

			product(q2, q1, normal);
			normal_mag = sqrt(iproduct(normal, normal));

			faces[nface].volume = fabs(iproduct(normal, dr)) / 6.0;
			faces[nface].area = normal_mag / 2.0;
			for (d = 0; d < DIMS; d++)
				faces[nface].normal[d] = normal[d] / normal_mag;

			// center of mass
			for (d = 0; d < DIMS; d++) {
				faces[nface].com_pbc[d] = (mons[n1].pos_pbc[d] + mons[n2].pos_pbc[d] + mons[n3].pos_pbc[d]) / 3.0;
				faces[nface].com[d] = box(faces[nface].com_pbc[d], maxsize[d]);
			}
		}
	}

	/* Normal vector of the monomers */
	for(n1 = 0; n1 < num_beads; n1++) {

		// Initialize
		for (d = 0; d < DIMS; d++)
			normal[d] = 0.0;
		normal_mag = 0.0;

		//
		for (j = 1; j <= mons[n1].Blist[0][0]; j++) {
			f2 = mons[n1].flist[j];

			for (d = 0; d < DIMS; d++)
				normal[d] += faces[f2].normal[d];
		}

		for (d = 0; d < DIMS; d++)
			normal[d] /= mons[n1].Blist[0][0];

		normal_mag = sqrt(iproduct(normal, normal));
		for (d = 0; d < DIMS; d++)
			mons[n1].normal[d] = normal[d] / normal_mag;
		mons[n1].curve = normal_mag;
	}

	free(delta_V);
	free(delta_A);
	for (i = 0; i < numsphere; i++)
		free(com[i]);
	free(com);
}

void sphere_props(struct sphere_param *sphere_pm, struct sphere *spheres, struct monomer *mono, struct face *faces)
{
	extern double maxsize[DIMS];
	int i, j, k, d, n1, n2, n3, d1, d2;

	int Nsphere = sphere_pm->Nsphere;
	int Ntype[NTYPES];
	int Nbead[NTYPES];
	int nspring, nface;
	int start, type;
	int min;
	int it_num, rot_num;
	double sprng_len, avg_sprng_len;
	double temp, tempstretch;
	double g_tensor[DIMS * DIMS];	/* gyration tensor */
	double I_tensor[DIMS * DIMS];	/* inertia tensor */
	double dr[DIMS], lambda[DIMS], Im[DIMS];
	double q12[DIMS], q13[DIMS], normal[DIMS];
	double y_com;
	double Rx2, Ry2, Rz2;
	double costheta, sintheta;
	double omega[DIMS];
	double y[162];
	double Iv[DIMS * DIMS];		/* syfan 151115 */
	double gv[DIMS * DIMS];
	FILE *stream;

	for (i = 0; i < NTYPES; i++) {
		Ntype[i] = sphere_pm->Ntype[i];
		Nbead[i] = sphere_pm->N_per_sphere[i];
	}

	// Calculate the sphere center-of-mass
	for (i = 0; i < Nsphere; i++) {
		type = (i < Ntype[0] ? 0 : 1);
		start = (type == 0 ? i * Nbead[0] : Ntype[0] * Nbead[0] + (i - Ntype[0]) * Nbead[1]);

		// Zero variables

		spheres[i].disp2 = 0.0;
		spheres[i].dr2 = 0.0;
		for (d = 0; d < DIMS; d++)
			spheres[i].com[d] = 0.0;

		// Calculate COM, the accumulative displacement(disp2), and the displacement between two adjacnt steps(dr2)

		for (j = 0; j < Nbead[type]; j++)
			for (d = 0; d < DIMS; d++)
				spheres[i].com[d] += mono[start + j].pos_pbc[d];
		for (d = 0; d < DIMS; d++) {
			spheres[i].com[d] = spheres[i].com[d] / Nbead[type];
			double temp1 = spheres[i].com[d] - spheres[i].com0[d];
			spheres[i].disp2 += (temp1 * temp1);
			double temp2 = spheres[i].com[d] - spheres[i].comold[d];
			spheres[i].dr2 += (temp2 * temp2);
			spheres[i].comold[d] = spheres[i].com[d];
		}
		spheres[i].disp2 = sqrt(spheres[i].disp2);
		spheres[i].dr2 = sqrt(spheres[i].dr2);
	}

	/* calculate the sphere radius of gyration and stretch */
	for (i = 0; i < Nsphere; i++) {
		type = (i < Ntype[0] ? 0 : 1);
		start = (type == 0 ? i * Nbead[0] : Ntype[0] * Nbead[0] + (i - Ntype[0]) * Nbead[1]);

		for (d = 0; d < DIMS; d++)
			spheres[i].stretch[d] = 0.0;

		for (d1 = 0; d1 < DIMS; d1++)
			for (d2 = 0; d2 < DIMS; d2++)
				g_tensor[d1 + DIMS * d2] = 0.0;

		for (j = 0; j < Nbead[type]; j++) {
			for (d = 0; d < DIMS; d++)
				dr[d] = mono[start + j].pos_pbc[d] - spheres[i].com[d];

			for (d1 = 0; d1 < DIMS; d1++)
				for (d2 = 0; d2 < DIMS; d2++)
					g_tensor[d1 + DIMS * d2] += dr[d1] * dr[d2];

			for (k = 0; k < Nbead[type]; k++) {
				tempstretch = 0.0;

				for (d = 0; d < DIMS; d++) {
					temp = mono[start + k].pos_pbc[d] - mono[start + j].pos_pbc[d];
					tempstretch += temp * temp;

					if (tempstretch > spheres[i].stretch[d])
						spheres[i].stretch[d] = tempstretch;
				}
			}
		}

		for (d1 = 0; d1 < DIMS; d1++)
			for (d2 = 0; d2 < DIMS; d2++)
				g_tensor[d1 + DIMS * d2] /= Nbead[type];

		/* convert the gyration tensor to the inertia tensor */
		I_tensor[0] = g_tensor[4] + g_tensor[8];
		I_tensor[4] = g_tensor[0] + g_tensor[8];
		I_tensor[8] = g_tensor[0] + g_tensor[4];
		I_tensor[1] = I_tensor[3] = -g_tensor[1];
		I_tensor[2] = I_tensor[6] = -g_tensor[2];
		I_tensor[5] = I_tensor[7] = -g_tensor[5];

		/* the moments of gyration tensor are the eigenvalues of gyration tensor */
		//sym33eigen_values(g_tensor, lambda);
		jacobi_eigenvalue(3, g_tensor, 10, gv, lambda, &it_num, &rot_num);

		/* principal moments of inertia */
		//sym33eigen_values(I_tensor, Im);
		jacobi_eigenvalue(3, I_tensor, 10, Iv, Im, &it_num, &rot_num);

		//printf("it_num = %d, rot_num = %d\n", it_num, rot_num);
		//printf("Iv1 = (%f, %f, %f)\tIv2 = (%le, %le, %le)\tIv3 = (%le, %le, %le)\n", Iv[0], Iv[1], Iv[2], Iv[3], Iv[4], Iv[5], Iv[6], Iv[7], Iv[8]);
		//printf("Im = (%f, %f, %f)\n", Im[0], Im[1], Im[2]);

		costheta = fabs(Im[2]);
		if (costheta > 1.)
			costheta = 1.;
		if (costheta < -1.)
			costheta = -1.;
		spheres[i].theta = acos(costheta) * 180. / M_PI;
		for (d = 0; d < DIMS; d++)
			spheres[i].eigenvector3[d] = Iv[2 + DIMS * d];

		Rx2 = lambda[0];
		Ry2 = lambda[1];
		Rz2 = lambda[2];

		spheres[i].Rg2 = Rx2 + Ry2 + Rz2;
		spheres[i].Rx2 = Rx2;
		spheres[i].Ry2 = Ry2;
		spheres[i].Rz2 = Rz2;
		// spheres[i].asphericity = Rz2 - (Rx2 + Ry2)/2.;
		spheres[i].acylindricity = Ry2 - Rx2;

		spheres[i].asphericity = ((Rx2 - Ry2) * (Rx2 - Ry2) + (Rx2 - Rz2) * (Rx2 - Rz2) + (Ry2 - Rz2) * (Ry2 - Rz2))
			/ (2. * (Rx2 + Ry2 + Rz2) * (Rx2 + Ry2 + Rz2));
		spheres[i].Ia = Im[0];
		spheres[i].Ib = Im[1];
		spheres[i].Ic = Im[2];
	}

	/* calculate average angular velocity */
	for (i = 0; i < Nsphere; i++) {
		type = (i < Ntype[0] ? 0 : 1);
		start = (type == 0 ? i * Nbead[0] : Ntype[0] * Nbead[0] + (i - Ntype[0]) * Nbead[1]);
		for (d = 0; d < DIMS; d++)
			spheres[i].omega[d] = 0.;
		for (j = 0; j < Nbead[type]; j++) {
			temp = 0.;
			for (d = 0; d < DIMS; d++) {
				dr[d] = mono[start + j].pos_pbc[d] - spheres[i].com[d];
				temp += dr[d] * dr[d];
			}
			product(dr, mono[start + j].vel, omega);
			for (d = 0; d < DIMS; d++)
				spheres[i].omega[d] += omega[d] / temp;
		}
		for (d = 0; d < DIMS; d++)
			spheres[i].omega[d] /= Nbead[type];
	}

	/* calculate the average spring length */
	for (i = 0; i < Nsphere; i++) {
		type = (i < Ntype[0] ? 0 : 1);
		start = (type == 0 ? i * Nbead[0] : Ntype[0] * Nbead[0] + (i - Ntype[0]) * Nbead[1]);
		nspring = 0;
		avg_sprng_len = 0.0;

		for (j = 0; j < Nbead[type]; j++)
			for (k = 1; k <= mono[start + j].blist[0][0]; k++) {
				n1 = start + j;
				n2 = mono[n1].blist[k][0];

				sprng_len = 0.0;
				for (d = 0; d < DIMS; d++) {
					temp = n_image(mono[n1].pos[d] - mono[n2].pos[d], maxsize[d]);
					sprng_len += temp * temp;
				}
				sprng_len = sqrt(sprng_len);

				avg_sprng_len += sprng_len;
				nspring++;
			}
		avg_sprng_len /= nspring;
		spheres[i].avg_sprng_len = avg_sprng_len;
	}

	/* calculate surface area and volume */
	for (i = 0; i < Nsphere; i++) {
		type = (i < Ntype[0] ? 0 : 1);

		start = (type == 0 ? i * sphere_pm->face_per_sphere[0] : Ntype[0] * sphere_pm->face_per_sphere[0] + (i - Ntype[0]) * sphere_pm->face_per_sphere[1]);

		spheres[i].area = 0.;
		spheres[i].volume = 0.;
		spheres[i].volume_dupin = 0.;

		/* loop over all faces */
		/* n1, n2, n3 are the three points of the triangle; q12, q13 are two sides; dr is vector from com to n1 */

		for (j = 0; j < sphere_pm->face_per_sphere[type]; j++) {
			double faceNormalMagnitude = 0.;
			double faceNormalMagnitude2 = 0.;

			nface = start + j;
			n1 = faces[nface].vertices[0];

			for (d = 0; d < DIMS; d++)
				dr[d] = mono[n1].pos_pbc[d] - spheres[i].com[d];

			n2 = faces[nface].vertices[1];

			for (d = 0; d < DIMS; d++)
				q12[d] = mono[n2].pos_pbc[d] - mono[n1].pos_pbc[d];

			n3 = faces[nface].vertices[2];

			for (d = 0; d < DIMS; d++)
				q13[d] = mono[n3].pos_pbc[d] - mono[n1].pos_pbc[d];

			product(q12, q13, normal);
			faceNormalMagnitude2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
			faceNormalMagnitude = sqrt(faceNormalMagnitude2);

			faces[nface].area = faceNormalMagnitude / 2.;

			spheres[i].area += faces[nface].area;

			spheres[i].volume += fabs(dr[0] * normal[0] + dr[1] * normal[1] + dr[2] * normal[2]) / 6.;

			y_com = (mono[n1].pos_pbc[1] + mono[n2].pos_pbc[1] + mono[n3].pos_pbc[1]) / 3.;
			/* make the y-component of the normal vector point outwards */
			if (y_com - spheres[i].com[1] >= 0.0)
				spheres[i].volume_dupin += fabs(normal[1] * y_com / 2.);
			else
				spheres[i].volume_dupin -= fabs(normal[1] * y_com / 2.);

			//          sphere[i].inertia[0] += (1./5.) * faces[nface].area * 
			//          (faceNormalMagnitude2 - normal[0]*normal[0]) * faceNormalMagnitude;
			//          
			//          sphere[i].inertia[1] += (1./5.) * faces[nface].area * 
			//          (- normal[0] * normal[1]) * faceNormalMagnitude;
			//
			//          sphere[i].inertia[2] += (1./5.) * faces[nface].area * 
			//          (- normal[0] * normal[2]) * faceNormalMagnitude;
			//
			//          sphere[i].inertia[4] += (1./5.) * faces[nface].area * 
			//          (faceNormalMagnitude2 - normal[1]*normal[1]) * faceNormalMagnitude;
			//
			//          sphere[i].inertia[5] += (1./5.) * faces[nface].area * 
			//          (- normal[1] * normal[2]) * faceNormalMagnitude;
			//
			//          sphere[i].inertia[8] += (1./5.) * faces[nface].area * 
			//          (faceNormalMagnitude2 - normal[2]*normal[2]) * faceNormalMagnitude;
		}
		//      sphere[i].inertia[3] = sphere[i].inertia[1];
		//      sphere[i].inertia[6] = sphere[i].inertia[2];
		//      sphere[i].inertia[7] = sphere[i].inertia[5];

	}
}
