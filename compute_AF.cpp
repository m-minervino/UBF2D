#include <vector>
#include <math.h>
#include "/usr/local/apps/tecplot2022r1/360ex_2022r1/include/TECIO.h"
#include "UBF2D.h"

using namespace std;

void compute_AF(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values,
				vector < vector < vector < vector < double > > > > &aux_values, vector < vector < vector < vector < double > > > > &aux_values_N,
				vector < vector < vector < vector < bool > > > > &flags, vector < vector < vector < vector < double > > > > &AF_elem,
				int t, int mD, int kD, int numVars, int compVarsN, int numAF, int numAF_ADD, int i_TOT, int64_t* kMax, bool NS, bool flux_scheme, int grad_scheme, bool ParForm, bool limit_integr,
				vector < bool > &near_ID, double rho_inf, double V_inf, double chord, int rho_id, int vx_id, int vz_id, int vgx_id, int vgz_id, vector < double > &X_POLE, vector < double > &Z_POLE,
				int U_FORM, bool FT_FORM, double x_CoR, double z_CoR, double ampl_x, double freq_x, double t0_x, double ampl_z, double freq_z, double t0_z,
				double ampl_a, double freq_a, double t0_a, bool MOV_GRD, vector < double > &time,
				double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				double a_1, double a_2, double a_3, double a_4, vector < vector < double > > &WLS, vector < vector < double > > &AERO) {
	//VOLUME integrals and surface integrals over OMEGA_SW boundary
	double VF3x=0;  //Explicit viscous force (volume itegral)
	double VF3z=0;
	for (int j = 0; j <numAF; j++) {
		for (int k = 0; k <kD; k++) {
			for (int m = mD; m <i_TOT-1-mD; m++) {
				if (limit_integr) {
					if ((j<20)||(j>31)) {	//Exclude those integrands without the mid-field property
						if (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]) {	//Also include the SW region (flags[][][0]) to correctly compute the modified compr. corr.
							AERO[j][t]=AERO[j][t]+AF_elem[m][k][j][t];
						}
					} else {
						AERO[j][t]=AERO[j][t]+AF_elem[m][k][j][t];
					}
				} else {
					AERO[j][t]=AERO[j][t]+AF_elem[m][k][j][t];
				}
			}
		}
	}
	//SURFACE integrals on VISCOUS region boundary
	for (int k = 0; k <kD; k++) {
		for (int m = mD; m <i_TOT-1-mD; m++) {
			if (flags[m][k][1][t]==true) {
				cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
								rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
				if (k != 0) {
					if (flags[m][k - 1][1][t] == false) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+22][t] = AERO[numAF+22][t] -
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+24][t] = AERO[numAF+24][t] + 
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+26][t]=AERO[numAF+26][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+22][t] = AERO[numAF+22][t] - 
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+24][t]=AERO[numAF+24][t] +
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) * (
									dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+26][t]=AERO[numAF+26][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*dz_1;
					}
				} else {
					if ((flags[i_TOT - 2 - m][k][1][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+22][t] = AERO[numAF+22][t] -
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+24][t] = AERO[numAF+24][t] + 
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+26][t]=AERO[numAF+26][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+22][t] = AERO[numAF+22][t] -
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+24][t]=AERO[numAF+24][t] +
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+26][t]=AERO[numAF+26][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*dz_1;
					}
				}
				if ((flags[m+1][k][1][t]==false)||((flags[m+1][k][1][t]==true)&&(m==i_TOT-2-mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+22][t] = AERO[numAF+22][t] - 
								rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 14][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 13][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 13][t]));
						if (m==i_TOT-2-mD) {
							AERO[numAF+28][t] = AERO[numAF+28][t] - 
									rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
									(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 14][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 13][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 13][t]));
						}
						//ONERA compressibility correction term
						AERO[numAF+24][t] = AERO[numAF+24][t] + 
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k][rho_id-1][t]+nodal_values[m + 1][k + 1][rho_id-1][t])/2)-rho_inf)*dz_2 +
								rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						if (m==i_TOT-2-mD) {
							AERO[numAF+30][t] = AERO[numAF+30][t] + 
									0.5*pow(V_inf,2)*(((nodal_values[m + 1][k][rho_id-1][t]+nodal_values[m + 1][k + 1][rho_id-1][t])/2)-rho_inf)*dz_2 +
									rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
									(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t] +
								rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+22][t] = AERO[numAF+22][t] - 
								rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k + 1][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k + 1][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
						if (m==i_TOT-2-mD) {
							AERO[numAF+28][t] = AERO[numAF+28][t] - 
									rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
									(dz_2 * ((aux_values_N[m + 1][k + 1][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
									dx_2 * ((aux_values_N[m + 1][k + 1][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
						}
						//ONERA compressibility correction term
						AERO[numAF+24][t]=AERO[numAF+24][t] +
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
								rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						if (m==i_TOT-2-mD) {
							AERO[numAF+30][t] = AERO[numAF+30][t] +
									0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
									rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
									(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
									dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t] +
								rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+26][t]=AERO[numAF+26][t]-((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])/2)*
															(((nodal_values[m+1][k][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_2+
															((nodal_values[m+1][k][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_2));
					//Formulation C - 2nd surface integral
					AERO[numAF+26][t]=AERO[numAF+26][t]+((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*dz_2;
				}
				if ((flags[m][k+1][1][t]==false)||((flags[m][k+1][1][t]==true)&&(k==kD-1))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+22][t] = AERO[numAF+22][t] - 
								rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 14][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 13][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+24][t] = AERO[numAF+24][t] + 
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k + 1][rho_id-1][t]+nodal_values[m][k + 1][rho_id-1][t])/2)-rho_inf)*dz_3 +
								rz_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t] +
								rz_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+22][t] = AERO[numAF+22][t] -
								rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * ((aux_values_N[m][k + 1][20][t] + aux_values_N[m + 1][k + 1][20][t]) / 2) +
								dx_3 * ((aux_values_N[m][k + 1][19][t] + aux_values_N[m + 1][k + 1][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+24][t]=AERO[numAF+24][t] +
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
								rz_3 * ((nodal_values[m+1][k+1][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m+1][k+1][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m+1][k+1][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t] +
								rz_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+26][t]=AERO[numAF+26][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])/2)*
															(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_3+
															((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_3));
					//Formulation C - 2nd surface integral
					AERO[numAF+26][t]=AERO[numAF+26][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*dz_3;
				}
				if ((flags[m-1][k][1][t]==false)||((flags[m-1][k][1][t]==true)&&(m==mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+22][t] = AERO[numAF+22][t] - 
								rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 14][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 13][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 13][t]));
						if (m==mD) {
							AERO[numAF+28][t] = AERO[numAF+28][t] - 
									rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
									(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 14][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 13][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 13][t]));
						}
						//ONERA compressibility correction term
						AERO[numAF+24][t] = AERO[numAF+24][t] + 
								0.5*pow(V_inf,2)*(((nodal_values[m][k + 1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
								rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						if (m==mD) {
							AERO[numAF+30][t] = AERO[numAF+30][t] + 
									0.5*pow(V_inf,2)*(((nodal_values[m][k + 1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
									rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
									(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t] +
								rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+22][t] = AERO[numAF+22][t] - 
								rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k][20][t] + aux_values_N[m][k + 1][20][t]) / 2) +
								dx_4 * ((aux_values_N[m][k][19][t] + aux_values_N[m][k + 1][19][t]) / 2));
						if (m==mD) {
							AERO[numAF+28][t] = AERO[numAF+28][t] - 
									rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
									(dz_4 * ((aux_values_N[m][k][20][t] + aux_values_N[m][k + 1][20][t]) / 2) +
									dx_4 * ((aux_values_N[m][k][19][t] + aux_values_N[m][k + 1][19][t]) / 2));
						}
						//ONERA compressibility correction term
						AERO[numAF+24][t]=AERO[numAF+24][t] +
								0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
								rz_4 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m][k][17][t]) / 2));
						if (m==mD) {
							AERO[numAF+30][t] = AERO[numAF+30][t] +
									0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
									rz_4 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
									(dz_4 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m][k][18][t]) / 2) +
									dx_4 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m][k][17][t]) / 2));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+26][t]=AERO[numAF+26][t] +
								rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m][k][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+26][t]=AERO[numAF+26][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])/2)*
															(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m][k][vx_id-1][t])/2)*dz_4+
															((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m][k][vz_id-1][t])/2)*(-dx_4));
					//Formulation C - 2nd surface integral
					AERO[numAF+26][t]=AERO[numAF+26][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t])/2)*dz_4;
				}
				if (k != 0) {
					if (flags[m][k - 1][1][t] == false) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+23][t] = AERO[numAF+23][t] + 
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+25][t]=AERO[numAF+25][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m + 1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+27][t]=AERO[numAF+27][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+23][t] = AERO[numAF+23][t] + 
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+25][t]=AERO[numAF+25][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+27][t]=AERO[numAF+27][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*(-dx_1);
					}
				} else {
					if ((flags[i_TOT - 2 - m][k][1][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+23][t] = AERO[numAF+23][t] +
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+25][t]=AERO[numAF+25][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m + 1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+27][t]=AERO[numAF+27][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+23][t] = AERO[numAF+23][t] +
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+25][t]=AERO[numAF+25][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+27][t]=AERO[numAF+27][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*(-dx_1);
					}
				}
				if ((flags[m+1][k][1][t]==false)||((flags[m+1][k][1][t]==true)&&(m==i_TOT-2-mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+23][t] = AERO[numAF+23][t] +
								rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 14][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 13][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 13][t]));
						if (m==i_TOT-2-mD) {
							AERO[numAF+29][t] = AERO[numAF+29][t] +
									rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
									(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 14][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 13][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 13][t]));
						}
						//ONERA compressibility correction term
						AERO[numAF+25][t]=AERO[numAF+25][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k][rho_id-1][t]+nodal_values[m + 1][k + 1][rho_id-1][t])/2)-rho_inf)*dx_2 -
								rx_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						if (m==i_TOT-2-mD) {
							AERO[numAF+31][t] = AERO[numAF+31][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m + 1][k][rho_id-1][t]+nodal_values[m + 1][k + 1][rho_id-1][t])/2)-rho_inf)*dx_2 -
									rx_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
									(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t] -
								rx_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+23][t] = AERO[numAF+23][t] +
								rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k + 1][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k + 1][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
						if (m==i_TOT-2-mD) {
							AERO[numAF+29][t] = AERO[numAF+29][t] +
									rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
									(dz_2 * ((aux_values_N[m + 1][k + 1][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
									dx_2 * ((aux_values_N[m + 1][k + 1][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
						}
						//ONERA compressibility correction term
						AERO[numAF+25][t]=AERO[numAF+25][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_2 -
								rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						if (m==i_TOT-2-mD) {
							AERO[numAF+31][t] = AERO[numAF+31][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_2 -
									rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
									(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
									dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t] -
								rx_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+27][t]=AERO[numAF+27][t]-((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t])/2)*
																(((nodal_values[m+1][k][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_2+
																((nodal_values[m+1][k][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_2));
					//Formulation C - 2nd surface integral
					AERO[numAF+27][t]=AERO[numAF+27][t]+((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*(-dx_2);
				}
				if ((flags[m][k+1][1][t]==false)||((flags[m][k+1][1][t]==true)&&(k==kD-1))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+23][t] = AERO[numAF+23][t] +
								rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 14][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 13][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+25][t]=AERO[numAF+25][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k + 1][rho_id-1][t]+nodal_values[m][k + 1][rho_id-1][t])/2)-rho_inf)*dx_3 -
								rx_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t] -
								rx_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+23][t] = AERO[numAF+23][t] +
								rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * ((aux_values_N[m][k + 1][20][t] + aux_values_N[m + 1][k + 1][20][t]) / 2) +
								dx_3 * ((aux_values_N[m][k + 1][19][t] + aux_values_N[m + 1][k + 1][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+25][t]=AERO[numAF+25][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
								rx_3 * ((nodal_values[m+1][k+1][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m+1][k+1][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m+1][k+1][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t] -
								rx_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+27][t]=AERO[numAF+27][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t])/2)*
																(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_3));
					//Formulation C - 2nd surface integral
					AERO[numAF+27][t]=AERO[numAF+27][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*(-dx_3);
				}
				if ((flags[m-1][k][1][t]==false)||((flags[m-1][k][1][t]==true)&&(m==mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+23][t] = AERO[numAF+23][t] +
								rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 14][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 13][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 13][t]));
						if (m==mD) {
						AERO[numAF+29][t] = AERO[numAF+29][t] +
								rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 14][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 13][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 13][t]));
						}
						//ONERA compressibility correction term
						AERO[numAF+25][t]=AERO[numAF+25][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m][k + 1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dx_4 -
								rx_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						if (m==mD) {
							AERO[numAF+31][t] = AERO[numAF+31][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k + 1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dx_4 -
									rx_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
									(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t] -
								rx_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+23][t] = AERO[numAF+23][t] +
								rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k][20][t] + aux_values_N[m][k + 1][20][t]) / 2) +
								dx_4 * ((aux_values_N[m][k][19][t] + aux_values_N[m][k + 1][19][t]) / 2));
						if (m==mD) {
							AERO[numAF+29][t] = AERO[numAF+29][t] +
									rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
									(dz_4 * ((aux_values_N[m][k][20][t] + aux_values_N[m][k + 1][20][t]) / 2) +
									dx_4 * ((aux_values_N[m][k][19][t] + aux_values_N[m][k + 1][19][t]) / 2));
						}
						//ONERA compressibility correction term
						AERO[numAF+25][t]=AERO[numAF+25][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dx_4 -
								rx_4 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m][k][17][t]) / 2));
						if (m==mD) {
							AERO[numAF+31][t] = AERO[numAF+31][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dx_4 -
									rx_4 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
									(dz_4 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m][k][18][t]) / 2) +
									dx_4 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m][k][17][t]) / 2));
						}
						//Formulation C - 3rd surface integral
						AERO[numAF+27][t]=AERO[numAF+27][t] -
								rx_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m][k][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+27][t]=AERO[numAF+27][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m][k][vx_id-1][t])/2)*dz_4+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m][k][vz_id-1][t])/2)*(-dx_4));
					//Formulation C - 2nd surface integral
					AERO[numAF+27][t]=AERO[numAF+27][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t])/2)*(-dx_4);
				}
			}
		}
	}
	//SURFACE integrals on WAVE region boundary
	for (int k = 0; k <kD; k++) {
		for (int m = mD; m <i_TOT-1-mD; m++) {
			if (flags[m][k][2][t]==true) {
				cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
								rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
				if (k != 0) {
					if (flags[m][k - 1][2][t] == false) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+38][t] = AERO[numAF+38][t] -
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+40][t] = AERO[numAF+40][t] + 
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+42][t]=AERO[numAF+42][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+38][t] = AERO[numAF+38][t] - 
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+40][t]=AERO[numAF+40][t] +
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) * (
									dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+42][t]=AERO[numAF+42][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*dz_1;
					}
				} else {
					if ((flags[i_TOT - 2 - m][k][2][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+38][t] = AERO[numAF+38][t] -
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+40][t] = AERO[numAF+40][t] + 
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+42][t]=AERO[numAF+42][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+38][t] = AERO[numAF+38][t] -
									rz_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+40][t]=AERO[numAF+40][t] +
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+42][t]=AERO[numAF+42][t] +
									rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*dz_1;
					}
				}
				if ((flags[m+1][k][2][t]==false)||((flags[m+1][k][2][t]==true)&&(m==i_TOT-2-mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+38][t] = AERO[numAF+38][t] - 
								rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 14][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 13][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+40][t] = AERO[numAF+40][t] + 
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k][rho_id-1][t]+nodal_values[m + 1][k + 1][rho_id-1][t])/2)-rho_inf)*dz_2 +
								rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t] +
								rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+38][t] = AERO[numAF+38][t] - 
								rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k + 1][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k + 1][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+40][t]=AERO[numAF+40][t] +
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
								rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t] +
								rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+42][t]=AERO[numAF+42][t]-((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])/2)*
															(((nodal_values[m+1][k][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_2+
															((nodal_values[m+1][k][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_2));
					//Formulation C - 2nd surface integral
					AERO[numAF+42][t]=AERO[numAF+42][t]+((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*dz_2;
				}
				if ((flags[m][k+1][2][t]==false)||((flags[m][k+1][2][t]==true)&&(k==kD-1))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+38][t] = AERO[numAF+38][t] - 
								rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 14][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 13][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+40][t] = AERO[numAF+40][t] + 
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k + 1][rho_id-1][t]+nodal_values[m][k + 1][rho_id-1][t])/2)-rho_inf)*dz_3 +
								rz_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t] +
								rz_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+38][t] = AERO[numAF+38][t] -
								rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * ((aux_values_N[m][k + 1][20][t] + aux_values_N[m + 1][k + 1][20][t]) / 2) +
								dx_3 * ((aux_values_N[m][k + 1][19][t] + aux_values_N[m + 1][k + 1][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+40][t]=AERO[numAF+40][t] +
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
								rz_3 * ((nodal_values[m+1][k+1][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m+1][k+1][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m+1][k+1][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t] +
								rz_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+42][t]=AERO[numAF+42][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])/2)*
															(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_3+
															((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_3));
					//Formulation C - 2nd surface integral
					AERO[numAF+42][t]=AERO[numAF+42][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*dz_3;
				}
				if ((flags[m-1][k][2][t]==false)||((flags[m-1][k][2][t]==true)&&(m==mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+38][t] = AERO[numAF+38][t] - 
								rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 14][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 13][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+40][t] = AERO[numAF+40][t] + 
								0.5*pow(V_inf,2)*(((nodal_values[m][k + 1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
								rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t] +
								rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+38][t] = AERO[numAF+38][t] - 
								rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k][20][t] + aux_values_N[m][k + 1][20][t]) / 2) +
								dx_4 * ((aux_values_N[m][k][19][t] + aux_values_N[m][k + 1][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+40][t]=AERO[numAF+40][t] +
								0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
								rz_4 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m][k][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+42][t]=AERO[numAF+42][t] +
								rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m][k][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+42][t]=AERO[numAF+42][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])/2)*
															(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m][k][vx_id-1][t])/2)*dz_4+
															((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m][k][vz_id-1][t])/2)*(-dx_4));
					//Formulation C - 2nd surface integral
					AERO[numAF+42][t]=AERO[numAF+42][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t])/2)*dz_4;
				}
				if (k != 0) {
					if (flags[m][k - 1][2][t] == false) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+39][t] = AERO[numAF+39][t] + 
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+41][t]=AERO[numAF+41][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m + 1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+43][t]=AERO[numAF+43][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+39][t] = AERO[numAF+39][t] + 
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+41][t]=AERO[numAF+41][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+43][t]=AERO[numAF+43][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*(-dx_1);
					}
				} else {
					if ((flags[i_TOT - 2 - m][k][2][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
						if (flux_scheme) {
							//Lamb vector circulation moment
							AERO[numAF+39][t] = AERO[numAF+39][t] +
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 14][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 14][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 13][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 13][t]));
							//ONERA compressibility correction term
							AERO[numAF+41][t]=AERO[numAF+41][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m + 1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
							//Formulation C - 3rd surface integral
							AERO[numAF+43][t]=AERO[numAF+43][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
									dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						} else {
							//Lamb vector circulation moment
							AERO[numAF+39][t] = AERO[numAF+39][t] +
									rx_1 * ((nodal_values[m + 1][k][rho_id-1][t] + nodal_values[m][k][rho_id-1][t]) / 2) * (
									dz_1 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_1 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m][k][19][t]) / 2));
							//ONERA compressibility correction term
							AERO[numAF+41][t]=AERO[numAF+41][t] -
									0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dx_1 -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
							//Formulation C - 3rd surface integral
							AERO[numAF+43][t]=AERO[numAF+43][t] -
									rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
									(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
									dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						}
						//Formulation C - 1st surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						//Formulation C - 2nd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*(-dx_1);
					}
				}
				if ((flags[m+1][k][2][t]==false)||((flags[m+1][k][2][t]==true)&&(m==i_TOT-2-mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+39][t] = AERO[numAF+39][t] +
								rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 14][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 13][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+41][t]=AERO[numAF+41][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k][rho_id-1][t]+nodal_values[m + 1][k + 1][rho_id-1][t])/2)-rho_inf)*dx_2 -
								rx_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t] -
								rx_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+39][t] = AERO[numAF+39][t] +
								rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k + 1][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k + 1][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+41][t]=AERO[numAF+41][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_2 -
								rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t] -
								rx_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
								(dz_2 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
								dx_2 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+43][t]=AERO[numAF+43][t]-((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t])/2)*
																(((nodal_values[m+1][k][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_2+
																((nodal_values[m+1][k][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_2));
					//Formulation C - 2nd surface integral
					AERO[numAF+43][t]=AERO[numAF+43][t]+((nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*(-dx_2);
				}
				if ((flags[m][k+1][2][t]==false)||((flags[m][k+1][2][t]==true)&&(k==kD-1))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+39][t] = AERO[numAF+39][t] +
								rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 14][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 13][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+41][t]=AERO[numAF+41][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m + 1][k + 1][rho_id-1][t]+nodal_values[m][k + 1][rho_id-1][t])/2)-rho_inf)*dx_3 -
								rx_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t] -
								rx_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+39][t] = AERO[numAF+39][t] +
								rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
								(dz_3 * ((aux_values_N[m][k + 1][20][t] + aux_values_N[m + 1][k + 1][20][t]) / 2) +
								dx_3 * ((aux_values_N[m][k + 1][19][t] + aux_values_N[m + 1][k + 1][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+41][t]=AERO[numAF+41][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
								rx_3 * ((nodal_values[m+1][k+1][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m+1][k+1][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m+1][k+1][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t] -
								rx_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
								(dz_3 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
								dx_3 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+43][t]=AERO[numAF+43][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t])/2)*
																(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_3));
					//Formulation C - 2nd surface integral
					AERO[numAF+43][t]=AERO[numAF+43][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*(-dx_3);
				}
				if ((flags[m-1][k][2][t]==false)||((flags[m-1][k][2][t]==true)&&(m==mD))) {
					if (flux_scheme) {
						//Lamb vector circulation moment
						AERO[numAF+39][t] = AERO[numAF+39][t] +
								rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 14][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 14][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 13][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 13][t]));
						//ONERA compressibility correction term
						AERO[numAF+41][t]=AERO[numAF+41][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m][k + 1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dx_4 -
								rx_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t] -
								rx_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
								dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						//Lamb vector circulation moment
						AERO[numAF+39][t] = AERO[numAF+39][t] +
								rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k][20][t] + aux_values_N[m][k + 1][20][t]) / 2) +
								dx_4 * ((aux_values_N[m][k][19][t] + aux_values_N[m][k + 1][19][t]) / 2));
						//ONERA compressibility correction term
						AERO[numAF+41][t]=AERO[numAF+41][t] -
								0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dx_4 -
								rx_4 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m][k][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+43][t]=AERO[numAF+43][t] -
								rx_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
								(dz_4 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m][k][18][t]) / 2) +
								dx_4 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m][k][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+43][t]=AERO[numAF+43][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m][k][vx_id-1][t])/2)*dz_4+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m][k][vz_id-1][t])/2)*(-dx_4));
					//Formulation C - 2nd surface integral
					AERO[numAF+43][t]=AERO[numAF+43][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t])/2)*(-dx_4);
				}
			}
		}
	}
	//SURFACE integrals over FAR-FIELD control volume boundary
	for (int k = 0; k <kD; k++) {
		if (k==kD-1) {
			for (int m = mD; m <i_TOT-1-mD; m++) {
				cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
							rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
				if (m==mD) {
					if (flux_scheme) {
						if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
							//Lamb vector circulation moment
							AERO[numAF][t]=AERO[numAF][t]-rz_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																(dz_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+14][t])+
																dx_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																	(dz_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+14][t])+
																	dx_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF][t]=AERO[numAF][t]-rz_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
																(dz_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+14][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+14][t])+
																dx_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+13][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
																	(dz_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+14][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+14][t])+
																	dx_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+13][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+13][t]));
							//Circulation moment of the viscous stress tensor divergence
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_3*
																(dz_3*(a_3 * aux_values[m][k+1][8][t] + (1 - a_3) * aux_values[m][k][8][t])+
																dx_3*(a_3 * aux_values[m][k+1][7][t] + (1 - a_3) * aux_values[m][k][7][t]));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_3*
																	(dz_3*(a_3 * aux_values[m][k+1][8][t] + (1 - a_3) * aux_values[m][k][8][t])+
																	dx_3*(a_3 * aux_values[m][k+1][7][t] + (1 - a_3) * aux_values[m][k][7][t]));
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_4*
																(dz_4*(a_4 * aux_values[m-1][k][8][t] + (1 - a_4) * aux_values[m][k][8][t])+
																dx_4*(a_4 * aux_values[m-1][k][7][t] + (1 - a_4) * aux_values[m][k][7][t]));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_4*
																	(dz_4*(a_4 * aux_values[m-1][k][8][t] + (1 - a_4) * aux_values[m][k][8][t])+
																	dx_4*(a_4 * aux_values[m-1][k][7][t] + (1 - a_4) * aux_values[m][k][7][t]));
							//Viscous stress tensor flux
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_3*(a_3 * aux_values[m][k+1][3][t] + (1 - a_3) * aux_values[m][k][3][t])-
																	dx_3*(a_3 * aux_values[m][k+1][4][t] + (1 - a_3) * aux_values[m][k][4][t]));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_3*(a_3 * aux_values[m][k+1][4][t] + (1 - a_3) * aux_values[m][k][4][t])-
																	dx_3*(a_3 * aux_values[m][k+1][6][t] + (1 - a_3) * aux_values[m][k][6][t]));
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_4*(a_4 * aux_values[m-1][k][3][t] + (1 - a_4) * aux_values[m][k][3][t])-
																	dx_4*(a_4 * aux_values[m-1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_4*(a_4 * aux_values[m-1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t])-
																	dx_4*(a_4 * aux_values[m-1][k][6][t] + (1 - a_4) * aux_values[m][k][6][t]));
						}
						//ONERA compressibility correction term
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dz_4 +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dx_4 -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Compressibility correction (Liu formulation - Surface integral term)
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
														(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
														(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
							//Lamb vector circulation moment
							AERO[numAF][t]=AERO[numAF][t]-rz_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																(dz_3*((aux_values_N[m][k+1][20][t]+aux_values_N[m+1][k+1][20][t])/2)+
																dx_3*((aux_values_N[m][k+1][19][t]+aux_values_N[m+1][k+1][19][t])/2));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																	(dz_3*((aux_values_N[m][k+1][20][t]+aux_values_N[m+1][k+1][20][t])/2)+
																	dx_3*((aux_values_N[m][k+1][19][t]+aux_values_N[m+1][k+1][19][t])/2));
							AERO[numAF][t]=AERO[numAF][t]-rz_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
																(dz_4*((aux_values_N[m][k][20][t]+aux_values_N[m][k+1][20][t])/2)+
																dx_4*((aux_values_N[m][k][19][t]+aux_values_N[m][k+1][19][t])/2));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
																	(dz_4*((aux_values_N[m][k][20][t]+aux_values_N[m][k+1][20][t])/2)+
																	dx_4*((aux_values_N[m][k][19][t]+aux_values_N[m][k+1][19][t])/2));
							//Circulation moment of the viscous stress tensor divergence
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_3*
																(dz_3*((aux_values_N[m+1][k+1][13][t]+aux_values_N[m][k+1][13][t])/2)+
																dx_3*((aux_values_N[m+1][k+1][12][t]+aux_values_N[m][k+1][12][t])/2));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_3*
																(dz_3*((aux_values_N[m+1][k+1][13][t]+aux_values_N[m][k+1][13][t])/2)+
																dx_3*((aux_values_N[m+1][k+1][12][t]+aux_values_N[m][k+1][12][t])/2));
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_4*
																(dz_4*((aux_values_N[m][k+1][13][t]+aux_values_N[m][k][13][t])/2)+
																dx_4*((aux_values_N[m][k+1][12][t]+aux_values_N[m][k][12][t])/2));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_4*
																(dz_4*((aux_values_N[m][k+1][13][t]+aux_values_N[m][k][13][t])/2)+
																dx_4*((aux_values_N[m][k+1][12][t]+aux_values_N[m][k][12][t])/2));
							//Viscous stress tensor flux
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_3*((aux_values_N[m+1][k+1][9][t]+aux_values_N[m][k+1][9][t])/2)-
																	dx_3*((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_3*((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2)-
																	dx_3*((aux_values_N[m+1][k+1][11][t]+aux_values_N[m][k+1][11][t])/2));
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_4*((aux_values_N[m][k+1][9][t]+aux_values_N[m][k][9][t])/2)-
																	dx_4*((aux_values_N[m][k+1][10][t]+aux_values_N[m][k][10][t])/2));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_4*((aux_values_N[m][k+1][10][t]+aux_values_N[m][k][10][t])/2)-
																	dx_4*((aux_values_N[m][k+1][11][t]+aux_values_N[m][k][11][t])/2));
						}
						//ONERA compressibility correction term
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
														rz_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dz_4 +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
														rx_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dx_4 -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						//Compressibility correction (Liu formulation - Surface integral term)
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
														(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
														(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+16][t]=AERO[numAF+16][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					AERO[numAF+16][t]=AERO[numAF+16][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_4+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_4));
					AERO[numAF+17][t]=AERO[numAF+17][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					AERO[numAF+17][t]=AERO[numAF+17][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t])/2)*
																(((nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_4+
																((nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_4));
					//Formulation C - 2nd surface integral
					AERO[numAF+18][t]=AERO[numAF+18][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*dz_3;
					AERO[numAF+18][t]=AERO[numAF+18][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*dz_4;
					AERO[numAF+19][t]=AERO[numAF+19][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*(-dx_3);
					AERO[numAF+19][t]=AERO[numAF+19][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
																nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*(-dx_4);
					//Alternative FT formulation: unsteady term on external control volume boundary
					if (FT_FORM) {
						if (U_FORM==0) {			//Use inertial time derivative of absolute fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
															(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
															dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
															(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
															dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
						} else if (U_FORM==1) {		//Use non-inertial time derivative of relative fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][26][t] + aux_values_N[m + 1][k + 1][26][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][25][t] + aux_values_N[m + 1][k + 1][25][t]) / 2));
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
															(dz_4 * ((aux_values_N[m][k][26][t] + aux_values_N[m][k + 1][26][t]) / 2) +
															dx_4 * ((aux_values_N[m][k][25][t] + aux_values_N[m][k + 1][25][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][26][t] + aux_values_N[m + 1][k + 1][26][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][25][t] + aux_values_N[m + 1][k + 1][25][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
															(dz_4 * ((aux_values_N[m][k][26][t] + aux_values_N[m][k + 1][26][t]) / 2) +
															dx_4 * ((aux_values_N[m][k][25][t] + aux_values_N[m][k + 1][25][t]) / 2));
						} else if (U_FORM==2) {		//Use non-inertial time derivative of relative fluid velocity
							//Formulation not applicable to the mixed inertial/non-inertial approach
						}
					}
					//Apparent force from relative motion term
					if (U_FORM==1) {
						AERO[numAF+36][t]=AERO[numAF+36][t] +
														rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
						AERO[numAF+36][t]=AERO[numAF+36][t] +
														rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
														(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
						AERO[numAF+37][t]=AERO[numAF+37][t] -
														rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
						AERO[numAF+37][t]=AERO[numAF+37][t] -
														rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
														(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
					}
				} else if (m==i_TOT-2-mD) {
					if (flux_scheme) {
						if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
							//Lamb vector circulation moment
							AERO[numAF][t]=AERO[numAF][t]-rz_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																(dz_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+14][t])+
																dx_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																	(dz_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+14][t])+
																	dx_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF][t]=AERO[numAF][t]-rz_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																(dz_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+14][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+14][t])+
																dx_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+13][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																(dz_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+14][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+14][t])+
																dx_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+13][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+13][t]));
							//Circulation moment of the viscous stress tensor divergence
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_3*
																(dz_3*(a_3 * aux_values[m][k+1][8][t] + (1 - a_3) * aux_values[m][k][8][t])+
																dx_3*(a_3 * aux_values[m][k+1][7][t] + (1 - a_3) * aux_values[m][k][7][t]));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_3*
																	(dz_3*(a_3 * aux_values[m][k+1][8][t] + (1 - a_3) * aux_values[m][k][8][t])+
																	dx_3*(a_3 * aux_values[m][k+1][7][t] + (1 - a_3) * aux_values[m][k][7][t]));
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_2*
																(dz_2*(a_2 * aux_values[m+1][k][8][t] + (1 - a_2) * aux_values[m][k][8][t])+
																dx_2*(a_2 * aux_values[m+1][k][7][t] + (1 - a_2) * aux_values[m][k][7][t]));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_2*
																	(dz_2*(a_2 * aux_values[m+1][k][8][t] + (1 - a_2) * aux_values[m][k][8][t])+
																	dx_2*(a_2 * aux_values[m+1][k][7][t] + (1 - a_2) * aux_values[m][k][7][t]));
							//Viscous stress tensor flux
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_3*(a_3 * aux_values[m][k+1][3][t] + (1 - a_3) * aux_values[m][k][3][t])-
																	dx_3*(a_3 * aux_values[m][k+1][4][t] + (1 - a_3) * aux_values[m][k][4][t]));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_3*(a_3 * aux_values[m][k+1][4][t] + (1 - a_3) * aux_values[m][k][4][t])-
																	dx_3*(a_3 * aux_values[m][k+1][6][t] + (1 - a_3) * aux_values[m][k][6][t]));
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_2*(a_2 * aux_values[m+1][k][3][t] + (1 - a_2) * aux_values[m][k][3][t])-
																	dx_2*(a_2 * aux_values[m+1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_2*(a_2 * aux_values[m+1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t])-
																	dx_2*(a_2 * aux_values[m+1][k][6][t] + (1 - a_2) * aux_values[m][k][6][t]));
						}
						//ONERA compressibility correction term
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_2 -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Compressibility correction (Liu formulation - Surface integral term)
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
														(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
														(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
							//Lamb vector circulation moment
							AERO[numAF][t]=AERO[numAF][t]-rz_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																(dz_3*((aux_values_N[m][k+1][20][t]+aux_values_N[m+1][k+1][20][t])/2)+
																dx_3*((aux_values_N[m][k+1][19][t]+aux_values_N[m+1][k+1][19][t])/2));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																	(dz_3*((aux_values_N[m][k+1][20][t]+aux_values_N[m+1][k+1][20][t])/2)+
																	dx_3*((aux_values_N[m][k+1][19][t]+aux_values_N[m+1][k+1][19][t])/2));
							AERO[numAF][t]=AERO[numAF][t]-rz_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																(dz_2*((aux_values_N[m+1][k+1][20][t]+aux_values_N[m+1][k][20][t])/2)+
																dx_2*((aux_values_N[m+1][k+1][19][t]+aux_values_N[m+1][k][19][t])/2));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_2*((aux_values_N[m+1][k+1][20][t]+aux_values_N[m+1][k][20][t])/2)+
																	dx_2*((aux_values_N[m+1][k+1][19][t]+aux_values_N[m+1][k][19][t])/2));
							//Circulation moment of the viscous stress tensor divergence
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_3*
																(dz_3*((aux_values_N[m+1][k+1][13][t]+aux_values_N[m][k+1][13][t])/2)+
																dx_3*((aux_values_N[m+1][k+1][12][t]+aux_values_N[m][k+1][12][t])/2));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_3*
																(dz_3*((aux_values_N[m+1][k+1][13][t]+aux_values_N[m][k+1][13][t])/2)+
																dx_3*((aux_values_N[m+1][k+1][12][t]+aux_values_N[m][k+1][12][t])/2));
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_2*
																(dz_2*((aux_values_N[m+1][k][13][t]+aux_values_N[m+1][k+1][13][t])/2)+
																dx_2*((aux_values_N[m+1][k][12][t]+aux_values_N[m+1][k+1][12][t])/2));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_2*
																(dz_2*((aux_values_N[m+1][k][13][t]+aux_values_N[m+1][k+1][13][t])/2)+
																dx_2*((aux_values_N[m+1][k][12][t]+aux_values_N[m+1][k+1][12][t])/2));
							//Viscous stress tensor flux
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_3*((aux_values_N[m+1][k+1][9][t]+aux_values_N[m][k+1][9][t])/2)-
																	dx_3*((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_3*((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2)-
																	dx_3*((aux_values_N[m+1][k+1][11][t]+aux_values_N[m][k+1][11][t])/2));
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_2*((aux_values_N[m+1][k][9][t]+aux_values_N[m+1][k+1][9][t])/2)-
																	dx_2*((aux_values_N[m+1][k][10][t]+aux_values_N[m+1][k+1][10][t])/2));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_2*((aux_values_N[m+1][k][10][t]+aux_values_N[m+1][k+1][10][t])/2)-
																	dx_2*((aux_values_N[m+1][k][11][t]+aux_values_N[m+1][k+1][11][t])/2));
						}
						//ONERA compressibility correction term
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
														rz_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
														rx_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_2 -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						//Compressibility correction (Liu formulation - Surface integral term)
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
														(dz_2 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
														dx_2 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
														(dz_2 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
														dx_2 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+16][t]=AERO[numAF+16][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					AERO[numAF+16][t]=AERO[numAF+16][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]+
																nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
																(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_2+
																((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_2));
					AERO[numAF+17][t]=AERO[numAF+17][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					AERO[numAF+17][t]=AERO[numAF+17][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]+
																nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
																(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_2+
																((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_2));
					//Formulation C - 2nd surface integral
					AERO[numAF+18][t]=AERO[numAF+18][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*dz_3;
					AERO[numAF+18][t]=AERO[numAF+18][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
																nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*dz_2;
					AERO[numAF+19][t]=AERO[numAF+19][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*(-dx_3);
					AERO[numAF+19][t]=AERO[numAF+19][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
																nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*(-dx_2);
					//Alternative FT formulation: unsteady term on external control volume boundary
					if (FT_FORM) {
						if (U_FORM==0) {			//Use inertial time derivative of absolute fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
							AERO[numAF+32][t]=AERO[numAF+32][t] -
														rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
														(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
														dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
														rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
														(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
														dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
						} else if (U_FORM==1) {		//Use non-inertial time derivative of relative fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][26][t] + aux_values_N[m + 1][k + 1][26][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][25][t] + aux_values_N[m + 1][k + 1][25][t]) / 2));
							AERO[numAF+32][t]=AERO[numAF+32][t] -
														rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
														(dz_2 * ((aux_values_N[m + 1][k + 1][26][t] + aux_values_N[m + 1][k][26][t]) / 2) +
														dx_2 * ((aux_values_N[m + 1][k + 1][25][t] + aux_values_N[m + 1][k][25][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][26][t] + aux_values_N[m + 1][k + 1][26][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][25][t] + aux_values_N[m + 1][k + 1][25][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
														rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
														(dz_2 * ((aux_values_N[m + 1][k + 1][26][t] + aux_values_N[m + 1][k][26][t]) / 2) +
														dx_2 * ((aux_values_N[m + 1][k + 1][25][t] + aux_values_N[m + 1][k][25][t]) / 2));
						} else if (U_FORM==2) {		//Use non-inertial time derivative of relative fluid velocity
							//Formulation not applicable to the mixed inertial/non-inertial approach
						}
					}
					//Apparent force from relative motion term
					if (U_FORM==1) {
						AERO[numAF+36][t]=AERO[numAF+36][t] +
														rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
						AERO[numAF+36][t]=AERO[numAF+36][t] +
													rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
													(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
													dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
						AERO[numAF+37][t]=AERO[numAF+37][t] -
														rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
						AERO[numAF+37][t]=AERO[numAF+37][t] -
													rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
													(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
													dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
					}
				} else {
					if (flux_scheme) {
						if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
							//Lamb vector circulation moment
							AERO[numAF][t]=AERO[numAF][t]-rz_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																(dz_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+14][t])+
																dx_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																	(dz_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+14][t])+
																	dx_3*(a_3 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+13][t]));
							//Circulation moment of the viscous stress tensor divergence
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_3*
																(dz_3*(a_3 * aux_values[m][k+1][8][t] + (1 - a_3) * aux_values[m][k][8][t])+
																dx_3*(a_3 * aux_values[m][k+1][7][t] + (1 - a_3) * aux_values[m][k][7][t]));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_3*
																	(dz_3*(a_3 * aux_values[m][k+1][8][t] + (1 - a_3) * aux_values[m][k][8][t])+
																	dx_3*(a_3 * aux_values[m][k+1][7][t] + (1 - a_3) * aux_values[m][k][7][t]));
							//Viscous stress tensor flux
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_3*(a_3 * aux_values[m][k+1][3][t] + (1 - a_3) * aux_values[m][k][3][t])-
																	dx_3*(a_3 * aux_values[m][k+1][4][t] + (1 - a_3) * aux_values[m][k][4][t]));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_3*(a_3 * aux_values[m][k+1][4][t] + (1 - a_3) * aux_values[m][k][4][t])-
																	dx_3*(a_3 * aux_values[m][k+1][6][t] + (1 - a_3) * aux_values[m][k][6][t]));
						}
						//ONERA compressibility correction term
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Compressibility correction (Liu formulation - Surface integral term)
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Formulation C - 3rd surface integral
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
					} else {
						if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
							//Lamb vector circulation moment
							AERO[numAF][t]=AERO[numAF][t]-rz_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																(dz_3*((aux_values_N[m][k+1][20][t]+aux_values_N[m+1][k+1][20][t])/2)+
																dx_3*((aux_values_N[m][k+1][19][t]+aux_values_N[m+1][k+1][19][t])/2));
							AERO[numAF+1][t]=AERO[numAF+1][t]+rx_3*((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*
																	(dz_3*((aux_values_N[m][k+1][20][t]+aux_values_N[m+1][k+1][20][t])/2)+
																	dx_3*((aux_values_N[m][k+1][19][t]+aux_values_N[m+1][k+1][19][t])/2));
							//Circulation moment of the viscous stress tensor divergence
							AERO[numAF+2][t]=AERO[numAF+2][t]+rz_3*
																(dz_3*((aux_values_N[m+1][k+1][13][t]+aux_values_N[m][k+1][13][t])/2)+
																dx_3*((aux_values_N[m+1][k+1][12][t]+aux_values_N[m][k+1][12][t])/2));
							AERO[numAF+3][t]=AERO[numAF+3][t]-rx_3*
																(dz_3*((aux_values_N[m+1][k+1][13][t]+aux_values_N[m][k+1][13][t])/2)+
																dx_3*((aux_values_N[m+1][k+1][12][t]+aux_values_N[m][k+1][12][t])/2));
							//Viscous stress tensor flux
							AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_3*((aux_values_N[m+1][k+1][9][t]+aux_values_N[m][k+1][9][t])/2)-
																	dx_3*((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2));
							AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_3*((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2)-
																	dx_3*((aux_values_N[m+1][k+1][11][t]+aux_values_N[m][k+1][11][t])/2));
						}
						//ONERA compressibility correction term
						AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
														rz_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_3 -
														rx_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						//Compressibility correction (Liu formulation - Surface integral term)
						AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_3 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
						//Formulation C - 3rd surface integral
						AERO[numAF+20][t]=AERO[numAF+20][t] +
														rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
						AERO[numAF+21][t]=AERO[numAF+21][t] -
														rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
					}
					//Formulation C - 1st surface integral
					AERO[numAF+16][t]=AERO[numAF+16][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					AERO[numAF+17][t]=AERO[numAF+17][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t])/2)*
																(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
																((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					//Formulation C - 2nd surface integral
					AERO[numAF+18][t]=AERO[numAF+18][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*dz_3;
					AERO[numAF+19][t]=AERO[numAF+19][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
																nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*(-dx_3);
					//Alternative FT formulation: unsteady term on external control volume boundary
					if (FT_FORM) {
						if (U_FORM==0) {			//Use inertial time derivative of absolute fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
						} else if (U_FORM==1) {		//Use non-inertial time derivative of relative fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t] -
															rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][26][t] + aux_values_N[m + 1][k + 1][26][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][25][t] + aux_values_N[m + 1][k + 1][25][t]) / 2));
							AERO[numAF+33][t]=AERO[numAF+33][t] +
															rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
															(dz_3 * ((aux_values_N[m][k + 1][26][t] + aux_values_N[m + 1][k + 1][26][t]) / 2) +
															dx_3 * ((aux_values_N[m][k + 1][25][t] + aux_values_N[m + 1][k + 1][25][t]) / 2));
						} else if (U_FORM==2) {		//Use non-inertial time derivative of relative fluid velocity
							//Formulation not applicable to the mixed inertial/non-inertial approach
						}
					}
					//Apparent force from relative motion term
					if (U_FORM==1) {
						AERO[numAF+36][t]=AERO[numAF+36][t] +
														rz_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
						AERO[numAF+37][t]=AERO[numAF+37][t] -
														rx_3 * ((nodal_values[m][k + 1][rho_id-1][t] + nodal_values[m + 1][k + 1][rho_id-1][t]) / 2) *
														(dz_3 * ((aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 2) +
														dx_3 * ((aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 2));
					}
				}
			}
		} else if (mD<i_TOT-1-mD) {
			int m=mD;
			cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
						rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
			if (flux_scheme) {
				if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
					//Lamb vector circulation moment
					AERO[numAF][t]=AERO[numAF][t]-rz_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
														(dz_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+14][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+14][t])+
														dx_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+13][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+13][t]));
					AERO[numAF+1][t]=AERO[numAF+1][t]+rx_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
															(dz_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+14][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+14][t])+
															dx_4*(a_4 * cell_values[m-1][k][numVars+compVarsN+13][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+13][t]));
					//Circulation moment of the viscous stress tensor divergence
					AERO[numAF+2][t]=AERO[numAF+2][t]+rz_4*
															(dz_4*(a_4 * aux_values[m-1][k][8][t] + (1 - a_4) * aux_values[m][k][8][t])+
															dx_4*(a_4 * aux_values[m-1][k][7][t] + (1 - a_4) * aux_values[m][k][7][t]));
					AERO[numAF+3][t]=AERO[numAF+3][t]-rx_4*
															(dz_4*(a_4 * aux_values[m-1][k][8][t] + (1 - a_4) * aux_values[m][k][8][t])+
															dx_4*(a_4 * aux_values[m-1][k][7][t] + (1 - a_4) * aux_values[m][k][7][t]));
					//Viscous stress tensor flux
					AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_4*(a_4 * aux_values[m-1][k][3][t] + (1 - a_4) * aux_values[m][k][3][t])-
																dx_4*(a_4 * aux_values[m-1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]));
					AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_4*(a_4 * aux_values[m-1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t])-
																dx_4*(a_4 * aux_values[m-1][k][6][t] + (1 - a_4) * aux_values[m][k][6][t]));
				}
				//ONERA compressibility correction term
				AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dz_4 +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
				AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dx_4 -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
				//Compressibility correction (Liu formulation - Surface integral term)
				AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
				AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
				//Formulation C - 3rd surface integral
				AERO[numAF+20][t]=AERO[numAF+20][t] +
												rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
												(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
												dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
				AERO[numAF+21][t]=AERO[numAF+21][t] -
												rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
												(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
												dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
			} else {
				if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
					//Lamb vector circulation moment
					AERO[numAF][t]=AERO[numAF][t]-rz_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
														(dz_4*((aux_values_N[m][k][20][t]+aux_values_N[m][k+1][20][t])/2)+
														dx_4*((aux_values_N[m][k][19][t]+aux_values_N[m][k+1][19][t])/2));
					AERO[numAF+1][t]=AERO[numAF+1][t]+rx_4*((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*
															(dz_4*((aux_values_N[m][k][20][t]+aux_values_N[m][k+1][20][t])/2)+
															dx_4*((aux_values_N[m][k][19][t]+aux_values_N[m][k+1][19][t])/2));
					//Circulation moment of the viscous stress tensor divergence
					AERO[numAF+2][t]=AERO[numAF+2][t]+rz_4*
														(dz_4*((aux_values_N[m][k+1][13][t]+aux_values_N[m][k][13][t])/2)+
															dx_4*((aux_values_N[m][k+1][12][t]+aux_values_N[m][k][12][t])/2));
					AERO[numAF+3][t]=AERO[numAF+3][t]-rx_4*
														(dz_4*((aux_values_N[m][k+1][13][t]+aux_values_N[m][k][13][t])/2)+
															dx_4*((aux_values_N[m][k+1][12][t]+aux_values_N[m][k][12][t])/2));
					//Viscous stress tensor flux
					AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_4*((aux_values_N[m][k+1][9][t]+aux_values_N[m][k][9][t])/2)-
																dx_4*((aux_values_N[m][k+1][10][t]+aux_values_N[m][k][10][t])/2));
					AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_4*((aux_values_N[m][k+1][10][t]+aux_values_N[m][k][10][t])/2)-
																dx_4*((aux_values_N[m][k+1][11][t]+aux_values_N[m][k][11][t])/2));
				}
				//ONERA compressibility correction term
				AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dz_4 +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
				AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dx_4 -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
				//Compressibility correction (Liu formulation - Surface integral term)
				AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
				AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
														dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k+1][17][t]) / 2));
				//Formulation C - 3rd surface integral
				AERO[numAF+20][t]=AERO[numAF+20][t] +
												rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
												(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
												dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
				AERO[numAF+21][t]=AERO[numAF+21][t] -
												rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
												(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
												dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
			}
			//Formulation C - 1st surface integral
			AERO[numAF+16][t]=AERO[numAF+16][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
														nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])/2)*
														(((nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_4+
														((nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_4));
			AERO[numAF+17][t]=AERO[numAF+17][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
														nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t])/2)*
														(((nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_4+
														((nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_4));
			//Formulation C - 2nd surface integral
			AERO[numAF+18][t]=AERO[numAF+18][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
														nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*dz_4;
			AERO[numAF+19][t]=AERO[numAF+19][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
														nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*(-dx_4);
			//Alternative FT formulation: unsteady term on external control volume boundary
			if (FT_FORM) {
				if (U_FORM==0) {			//Use inertial time derivative of absolute fluid velocity
					AERO[numAF+32][t]=AERO[numAF+32][t] -
													rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
													dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
					AERO[numAF+33][t]=AERO[numAF+33][t] +
													rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
													dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
				} else if (U_FORM==1) {		//Use non-inertial time derivative of relative fluid velocity
					AERO[numAF+32][t]=AERO[numAF+32][t] -
													rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][26][t] + aux_values_N[m][k + 1][26][t]) / 2) +
													dx_4 * ((aux_values_N[m][k][25][t] + aux_values_N[m][k + 1][25][t]) / 2));
					AERO[numAF+33][t]=AERO[numAF+33][t] +
													rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][26][t] + aux_values_N[m][k + 1][26][t]) / 2) +
													dx_4 * ((aux_values_N[m][k][25][t] + aux_values_N[m][k + 1][25][t]) / 2));
				} else if (U_FORM==2) {		//Use non-inertial time derivative of relative fluid velocity
					//Formulation not applicable to the mixed inertial/non-inertial approach
				}
			}
			//Apparent force from relative motion term
			if (U_FORM==1) {
				AERO[numAF+36][t]=AERO[numAF+36][t] +
												rz_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
												(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
												dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
				AERO[numAF+37][t]=AERO[numAF+37][t] -
												rx_4 * ((nodal_values[m][k][rho_id-1][t] + nodal_values[m][k + 1][rho_id-1][t]) / 2) *
												(dz_4 * ((aux_values_N[m][k][28][t] + aux_values_N[m][k + 1][28][t]) / 2) +
												dx_4 * ((aux_values_N[m][k][27][t] + aux_values_N[m][k + 1][27][t]) / 2));
			}
			m=i_TOT-2-mD;
			cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
						rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
			if (flux_scheme) {
				if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
					//Lamb vector circulation moment
					AERO[numAF][t]=AERO[numAF][t]-rz_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
														(dz_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+14][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+14][t])+
														dx_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+13][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+13][t]));
					AERO[numAF+1][t]=AERO[numAF+1][t]+rx_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
														(dz_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+14][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+14][t])+
														dx_2*(a_2 * cell_values[m+1][k][numVars+compVarsN+13][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+13][t]));
					//Circulation moment of the viscous stress tensor divergence
					AERO[numAF+2][t]=AERO[numAF+2][t]+rz_2*
															(dz_2*(a_2 * aux_values[m+1][k][8][t] + (1 - a_2) * aux_values[m][k][8][t])+
															dx_2*(a_2 * aux_values[m+1][k][7][t] + (1 - a_2) * aux_values[m][k][7][t]));
					AERO[numAF+3][t]=AERO[numAF+3][t]-rx_2*
															(dz_2*(a_2 * aux_values[m+1][k][8][t] + (1 - a_2) * aux_values[m][k][8][t])+
															dx_2*(a_2 * aux_values[m+1][k][7][t] + (1 - a_2) * aux_values[m][k][7][t]));
					//Viscous stress tensor flux
					AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_2*(a_2 * aux_values[m+1][k][3][t] + (1 - a_2) * aux_values[m][k][3][t])-
																dx_2*(a_2 * aux_values[m+1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]));
					AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_2*(a_2 * aux_values[m+1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t])-
																dx_2*(a_2 * aux_values[m+1][k][6][t] + (1 - a_2) * aux_values[m][k][6][t]));
				}
				//ONERA compressibility correction term
				AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
				AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_2 -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
				//Compressibility correction (Liu formulation - Surface integral term)
				AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
				AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
														dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
				//Formulation C - 3rd surface integral
				AERO[numAF+20][t]=AERO[numAF+20][t] +
												rz_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
												(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
												dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
				AERO[numAF+21][t]=AERO[numAF+21][t] -
												rx_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
												(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
												dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
			} else {
				if ((!limit_integr)||(limit_integr && (flags[m][k][0][t]||flags[m][k][1][t]||flags[m][k][2][t]))) {
					//Lamb vector circulation moment
					AERO[numAF][t]=AERO[numAF][t]-rz_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
														(dz_2*((aux_values_N[m+1][k+1][20][t]+aux_values_N[m+1][k][20][t])/2)+
														dx_2*((aux_values_N[m+1][k+1][19][t]+aux_values_N[m+1][k][19][t])/2));
					AERO[numAF+1][t]=AERO[numAF+1][t]+rx_2*((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
															(dz_2*((aux_values_N[m+1][k+1][20][t]+aux_values_N[m+1][k][20][t])/2)+
															dx_2*((aux_values_N[m+1][k+1][19][t]+aux_values_N[m+1][k][19][t])/2));
					//Circulation moment of the viscous stress tensor divergence
					AERO[numAF+2][t]=AERO[numAF+2][t]+rz_2*
														(dz_2*((aux_values_N[m+1][k][13][t]+aux_values_N[m+1][k+1][13][t])/2)+
															dx_2*((aux_values_N[m+1][k][12][t]+aux_values_N[m+1][k+1][12][t])/2));
					AERO[numAF+3][t]=AERO[numAF+3][t]-rx_2*
														(dz_2*((aux_values_N[m+1][k][13][t]+aux_values_N[m+1][k+1][13][t])/2)+
															dx_2*((aux_values_N[m+1][k][12][t]+aux_values_N[m+1][k+1][12][t])/2));
					//Viscous stress tensor flux
					AERO[numAF+4][t]=AERO[numAF+4][t]+(dz_2*((aux_values_N[m+1][k][9][t]+aux_values_N[m+1][k+1][9][t])/2)-
																dx_2*((aux_values_N[m+1][k][10][t]+aux_values_N[m+1][k+1][10][t])/2));
					AERO[numAF+5][t]=AERO[numAF+5][t]+(dz_2*((aux_values_N[m+1][k][10][t]+aux_values_N[m+1][k+1][10][t])/2)-
																dx_2*((aux_values_N[m+1][k][11][t]+aux_values_N[m+1][k+1][11][t])/2));
				}
				//ONERA compressibility correction term
				AERO[numAF+6][t]=AERO[numAF+6][t] + 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
				AERO[numAF+7][t]=AERO[numAF+7][t] - 0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dx_2 -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
				//Compressibility correction (Liu formulation - Surface integral term)
				AERO[numAF+14][t]=AERO[numAF+14][t] +
														rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
				AERO[numAF+15][t]=AERO[numAF+15][t] -
														rx_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k + 1][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
														dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
				//Formulation C - 3rd surface integral
				AERO[numAF+20][t]=AERO[numAF+20][t] +
												rz_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
												(dz_2 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
												dx_2 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
				AERO[numAF+21][t]=AERO[numAF+21][t] -
												rx_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
												(dz_2 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
												dx_2 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
			}
			//Formulation C - 1st surface integral
			AERO[numAF+16][t]=AERO[numAF+16][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]+
														nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
														(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_2+
														((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_2));
			AERO[numAF+17][t]=AERO[numAF+17][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]+
														nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
														(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_2+
														((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_2));
			//Formulation C - 2nd surface integral
			AERO[numAF+18][t]=AERO[numAF+18][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
														nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*dz_2;
			AERO[numAF+19][t]=AERO[numAF+19][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
														nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*(-dx_2);
			//Alternative FT formulation: unsteady term on external control volume boundary
			if (FT_FORM) {
				if (U_FORM==0) {			//Use inertial time derivative of absolute fluid velocity
					AERO[numAF+32][t]=AERO[numAF+32][t] -
												rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
												(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
												dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
					AERO[numAF+33][t]=AERO[numAF+33][t] +
												rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
												(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
												dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
				} else if (U_FORM==1) {		//Use non-inertial time derivative of relative fluid velocity
					AERO[numAF+32][t]=AERO[numAF+32][t] -
												rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
												(dz_2 * ((aux_values_N[m + 1][k + 1][26][t] + aux_values_N[m + 1][k][26][t]) / 2) +
												dx_2 * ((aux_values_N[m + 1][k + 1][25][t] + aux_values_N[m + 1][k][25][t]) / 2));
					AERO[numAF+33][t]=AERO[numAF+33][t] +
												rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
												(dz_2 * ((aux_values_N[m + 1][k + 1][26][t] + aux_values_N[m + 1][k][26][t]) / 2) +
												dx_2 * ((aux_values_N[m + 1][k + 1][25][t] + aux_values_N[m + 1][k][25][t]) / 2));
				} else if (U_FORM==2) {		//Use non-inertial time derivative of relative fluid velocity
					//Formulation not applicable to the mixed inertial/non-inertial approach
				}
			}
			//Apparent force from relative motion term
			if (U_FORM==1) {
				AERO[numAF+36][t]=AERO[numAF+36][t] +
											rz_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
											(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
											dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
				AERO[numAF+37][t]=AERO[numAF+37][t] -
											rx_2 * ((nodal_values[m + 1][k + 1][rho_id-1][t] + nodal_values[m + 1][k][rho_id-1][t]) / 2) *
											(dz_2 * ((aux_values_N[m + 1][k + 1][28][t] + aux_values_N[m + 1][k][28][t]) / 2) +
											dx_2 * ((aux_values_N[m + 1][k + 1][27][t] + aux_values_N[m + 1][k][27][t]) / 2));
			}
		}
		//SURFACE integrals over BODY boundary (only computed in case of eff. inviscid flow or moving-grid and inertial or mixed formulation)
		if ( (k==0) && ( ParForm || ((!NS)||(MOV_GRD && (U_FORM==0))) ) ) {		//The check on V_wall!=0 is done on the working velocity, for U_FORM=0 or U_FORM=1, while for U_FORM=2 it is done for the relative velocity.
			for (int m = mD; m <i_TOT-1-mD; m++) {
				if ((near_ID[m]) || (near_ID[m + 1])) {
					cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
								rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
					//Vortex force additional term on body surface
					if (U_FORM<2) {		//Use the working velocity (absolute for U_FORM=0, relative for U_FORM=1)
						AERO[numAF+8][t]=AERO[numAF+8][t]-((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2)*dz_1;
						AERO[numAF+9][t]=AERO[numAF+9][t]+((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2)*dx_1;
					} else {			//Use the relative velocity, as this is the velocity that appears in the F_body contribution of the mixed inertial/non-inertial formulation
						AERO[numAF+8][t]=AERO[numAF+8][t]-((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*((
								0.5*(pow(nodal_values[m][k][vx_id-1][t]-nodal_values[m][k][vgx_id-1][t],2)+pow(nodal_values[m][k][vz_id-1][t]-nodal_values[m][k][vgz_id-1][t],2))+
								0.5*(pow(nodal_values[m+1][k][vx_id-1][t]-nodal_values[m+1][k][vgx_id-1][t],2)+pow(nodal_values[m+1][k][vz_id-1][t]-nodal_values[m+1][k][vgz_id-1][t],2))
								) / 2)*dz_1;
						AERO[numAF+9][t]=AERO[numAF+9][t]+((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*((
								0.5*(pow(nodal_values[m][k][vx_id-1][t]-nodal_values[m][k][vgx_id-1][t],2)+pow(nodal_values[m][k][vz_id-1][t]-nodal_values[m][k][vgz_id-1][t],2))+
								0.5*(pow(nodal_values[m+1][k][vx_id-1][t]-nodal_values[m+1][k][vgx_id-1][t],2)+pow(nodal_values[m+1][k][vz_id-1][t]-nodal_values[m+1][k][vgz_id-1][t],2))
								) / 2)*dx_1;
					}
					if (flux_scheme) {
						//Compressibility correction additional term on body surface
						AERO[numAF+10][t]=AERO[numAF+10][t]-rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
																(dz_1 * (-a_1 * cell_values[m][k+1][numVars + compVarsN + 9][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
																dx_1 * (-a_1 * cell_values[m][k+1][numVars + compVarsN + 8][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						AERO[numAF+11][t]=AERO[numAF+11][t]+rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
																(dz_1 * (-a_1 * cell_values[m][k+1][numVars + compVarsN + 9][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
																dx_1 * (-a_1 * cell_values[m][k+1][numVars + compVarsN + 8][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
						//Outer vorticity contribution additional term on body surface
						if (!ParForm) {
							AERO[numAF+12][t]=AERO[numAF+12][t]+rz_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*(-a_1 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+14][t])+
																	dx_1*(-a_1 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+13][t]));
							AERO[numAF+13][t]=AERO[numAF+13][t]-rx_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*(-a_1 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+14][t])+
																	dx_1*(-a_1 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+13][t]));
						} else {
							AERO[numAF+12][t]=AERO[numAF+12][t]+rz_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((-a_1 * cell_values[m][k+1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t])*
																		(-(nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2))+
																	dx_1*((-a_1 * cell_values[m][k+1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t])*
																		((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)));
							AERO[numAF+13][t]=AERO[numAF+13][t]-rx_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((-a_1 * cell_values[m][k+1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t])*
																		(-(nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2))+
																	dx_1*((-a_1 * cell_values[m][k+1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t])*
																		((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)));
						}
					} else {
						//Compressibility correction additional term on body surface
						AERO[numAF+10][t]=AERO[numAF+10][t]-rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
																(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
																dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
						AERO[numAF+11][t]=AERO[numAF+11][t]+rx_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) *
																(dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
																dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
						//Outer vorticity contribution additional term on body surface
						if (!ParForm) {
							AERO[numAF+12][t]=AERO[numAF+12][t]+rz_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((aux_values_N[m][k][20][t]+aux_values_N[m+1][k][20][t])/2)+
																	dx_1*((aux_values_N[m][k][19][t]+aux_values_N[m+1][k][19][t])/2));
							AERO[numAF+13][t]=AERO[numAF+13][t]-rx_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((aux_values_N[m][k][20][t]+aux_values_N[m+1][k][20][t])/2)+
																	dx_1*((aux_values_N[m][k][19][t]+aux_values_N[m+1][k][19][t])/2));
						} else {
							AERO[numAF+12][t]=AERO[numAF+12][t]+rz_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((-(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vx_id-1][t]-
																			(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vx_id-1][t])/2)+
																	dx_1*(((aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vz_id-1][t]+
																			(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vz_id-1][t])/2));
							AERO[numAF+12][t]=AERO[numAF+12][t]-rx_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((-(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vx_id-1][t]-
																			(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vx_id-1][t])/2)+
																	dx_1*(((aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vz_id-1][t]+
																			(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vz_id-1][t])/2));
						}
					}
				}
			}
		}
		//SURFACE integrals over BODY boundary (always computed)
		if (k==0) {
			for (int m = mD; m <i_TOT-1-mD; m++) {
				if ((near_ID[m]) || (near_ID[m + 1])) {
					cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
								rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
					//Volume formulation of ONERA compressibility correction
					AERO[23][t]=AERO[23][t]-(1/rho_inf)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1;
					AERO[25][t]=AERO[25][t]-(1/rho_inf)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*(-dx_1);
					//Unsteady body surface term
					if (!FT_FORM) {
						if (U_FORM==0) {			//Use inertial time derivative of absolute fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t]+rz_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t])/2)+
																	dx_1*((aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t])/2));
							AERO[numAF+33][t]=AERO[numAF+33][t]-rx_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t])/2)+
																	dx_1*((aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t])/2));
						} else if (U_FORM==1) {			//Use non-inertial time derivative of relative fluid velocity
							AERO[numAF+32][t]=AERO[numAF+32][t]+rz_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((aux_values_N[m][k][26][t]+aux_values_N[m+1][k][26][t])/2)+
																	dx_1*((aux_values_N[m][k][25][t]+aux_values_N[m+1][k][25][t])/2));
							AERO[numAF+33][t]=AERO[numAF+33][t]-rx_1*((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
																	(dz_1*((aux_values_N[m][k][26][t]+aux_values_N[m+1][k][26][t])/2)+
																	dx_1*((aux_values_N[m][k][25][t]+aux_values_N[m+1][k][25][t])/2));
						} else if (U_FORM==2) {		//Use non-inertial time derivative of relative fluid velocity, plus moving reference terms
							double X_prime_x = (nodal_values[m][k][0][t]+nodal_values[m+1][k][0][t])/2 - x_CoR - ampl_x * sin(freq_x*(time[t]-t0_x));
							double X_prime_z = (nodal_values[m][k][2][t]+nodal_values[m+1][k][2][t])/2 - z_CoR - ampl_z * sin(freq_z*(time[t]-t0_z));
							double vx=NAN, vz=NAN;
							vx = (((aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t])/2) - freq_x*freq_x*ampl_x*sin(freq_x*(time[t]-t0_x)) -
								freq_a*freq_a*ampl_a*sin(freq_a*(time[t]-t0_a))*X_prime_z - pow(freq_a*ampl_a*cos(freq_a*(time[t]-t0_a)),2)*X_prime_x)*
								((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2);
							vz = (((aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t])/2) - freq_z*freq_z*ampl_z*sin(freq_z*(time[t]-t0_z)) +
								freq_a*freq_a*ampl_a*sin(freq_a*(time[t]-t0_a))*X_prime_x - pow(freq_a*ampl_a*cos(freq_a*(time[t]-t0_a)),2)*X_prime_z)*
								((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2);
							AERO[numAF+32][t]=AERO[numAF+32][t]+rz_1*(dz_1*vz + dx_1*vx);
							AERO[numAF+33][t]=AERO[numAF+33][t]-rx_1*(dz_1*vz + dx_1*vx);
						}
					}
					//Unsteady body surface term (Formulation C)
					if (MOV_GRD && (U_FORM!=1)) {	//Not computed in cases it should be exactely zero (to avoid inaccuracies due to non-zero relative velocity at trailing edge points)
						AERO[numAF+34][t]=AERO[numAF+34][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
																	(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																	((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
						AERO[numAF+35][t]=AERO[numAF+35][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
																	nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
																	(((nodal_values[m][k][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_1+
																	((nodal_values[m][k][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_1));
					}
				}
			}
		}
	}
	for (int var=0; var<numAF_ADD; var++) {		//SURFACE integrals non-dimensionalization
		if (mD<floor((i_TOT-1)/2)+1) {
			AERO[numAF+var][t] = AERO[numAF+var][t] /(0.5*rho_inf*pow(V_inf,2)*chord);
		}
	}
}