#include <vector>
#include "/opt/tecplot/360ex_2025r2/include/TECIO.h"

using namespace std;

void updateParasiteNew(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector <vector < double > > > > &cell_values,
					vector < double > &X_POLE, vector < double > &Z_POLE, int rho_id, int vx_id, int vz_id, vector < vector < vector <vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < bool > > > > &flags,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, bool MOV_GRD, bool flux_scheme, bool smart_wds,
					vector < bool > &near_ID, double CdP_tsh, double rho_inf, double V_inf, double chord,
					double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem) {
	//******** WARNING *********
	//This function computes gradients (and divergence) using the G-G nodes-based method only.
	//Other gradient schemes are to be included and implemented, as done in the rest of the code.
	//**************************
	double CdeP_vrt_1=0;
	double CdeP_vrt_2 = 0;
	//X-component
	if (flux_scheme) {
		//First term
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			CdeP_vrt_1 = 0;
		} else if (k != 0) {
			CdeP_vrt_1 = -((a_1 * cell_values[m][k - 1][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_1 *
							(((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 +
									(a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_2 *
							(((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 +
									(a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_3 *
							(((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 +
									(a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_4 *
							(((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				CdeP_vrt_1 = -((-a_1 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_1 *
							(((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 +
									(a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_2 *
							(((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 +
									(a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_3 *
							(((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 +
									(a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_4 *
							(((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
			} else {
				CdeP_vrt_1 = -((a_1 * cell_values[i_TOT - 2 - m][k][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_1 *
							(((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 +
									(a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_2 *
							(((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 +
									(a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_3 *
							(((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 +
									(a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t]) * dz_4 *
							(((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
			}
		}
		//Second term
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			CdeP_vrt_2 = 0;
		} else if (k != 0) {
			CdeP_vrt_2 = (
				(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * cell_values[m][k - 1][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) *
											(nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] + nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t]) * 0.5 +
				(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) *
											(nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] + nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]) * 0.5 +
				(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) *
											(nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] + nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]) * 0.5 +
				(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t])*
											(nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] + nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]) * 0.5);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				CdeP_vrt_2 = (
					(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] + nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t]) * 0.5 +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] + nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]) * 0.5 +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] + nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]) * 0.5 +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t])*
												(nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] + nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]) * 0.5);
			} else {
				CdeP_vrt_2 = (
					(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * cell_values[i_TOT - 2 - m][k][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] + nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t]) * 0.5 +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] + nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]) * 0.5 +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] + nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]) * 0.5 +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t])*
												(nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] + nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]) * 0.5);
			}
		}
	} else {
		//First term
		CdeP_vrt_1 = -( dz_1 * ((aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t]) *
								((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
								(nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
								(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t]) *
								((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
								(nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 +
						dz_2 * ((aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t]) *
								((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
								(nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
								(aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t]) *
								((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
								(nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 +
						dz_3 * ((aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t]) *
								((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
								(nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
								(aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t]) *
								((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
								(nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 +
						dz_4 * ((aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t]) *
								((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
								(nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
								(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t]) *
								((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
								(nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
		//Second term
		CdeP_vrt_2 = (
			(rx_1 * dz_1 - rz_1 * dx_1) * (nodal_values[m][k][rho_id-1][t]*(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vz_id-1][t] + 
											nodal_values[m+1][k][rho_id-1][t]*(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vz_id-1][t]) * 0.5 +
			(rx_2 * dz_2 - rz_2 * dx_2) * (nodal_values[m+1][k][rho_id-1][t]*(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vz_id-1][t] + 
											nodal_values[m+1][k+1][rho_id-1][t]*(aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t])*nodal_values[m+1][k+1][vz_id-1][t]) * 0.5 +
			(rx_3 * dz_3 - rz_3 * dx_3) * (nodal_values[m+1][k+1][rho_id-1][t]*(aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t])*nodal_values[m+1][k+1][vz_id-1][t] +
											nodal_values[m][k+1][rho_id-1][t]*(aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t])*nodal_values[m][k+1][vz_id-1][t]) * 0.5 +
			(rx_4 * dz_4 - rz_4 * dx_4) * (nodal_values[m][k+1][rho_id-1][t]*(aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t])*nodal_values[m][k+1][vz_id-1][t] +
											nodal_values[m][k][rho_id-1][t]*(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vz_id-1][t]) * 0.5);
	}
	//Attribute X-component first and second terms to different regions (WAVE, VISCOUS or SPURIOUS) based on total x-component value for the computed cell
	if (smart_wds) {
		if ((((CdeP_vrt_1+ CdeP_vrt_2) > CdP_tsh* 0.5 * rho_inf * V_inf * V_inf * chord) || ((CdeP_vrt_1 + CdeP_vrt_2) < -CdP_tsh * 0.5 * rho_inf * V_inf * V_inf * chord)) && flags[m][k][2][t]) {
			AF_elem[m][k][14][t] = CdeP_vrt_1;
			AF_elem[m][k][17][t] = CdeP_vrt_2;
		} else if (flags[m][k][1][t]) {
			AF_elem[m][k][15][t] = CdeP_vrt_1;
			AF_elem[m][k][18][t] = CdeP_vrt_2;
		} else {
			AF_elem[m][k][33][t] = CdeP_vrt_1;
			AF_elem[m][k][34][t] = CdeP_vrt_2;
			flags[m][k][2][t] = false;
		}
	} else {
		if (flags[m][k][2][t]) {
			AF_elem[m][k][14][t] = CdeP_vrt_1;
			AF_elem[m][k][17][t] = CdeP_vrt_2;
		} else if (flags[m][k][1][t]) {
			AF_elem[m][k][15][t] = CdeP_vrt_1;
			AF_elem[m][k][18][t] = CdeP_vrt_2;
		} else {
			AF_elem[m][k][33][t] = CdeP_vrt_1;
			AF_elem[m][k][34][t] = CdeP_vrt_2;
		}
	}
	//Z-component
	//---First term
	if (flux_scheme) {
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			AF_elem[m][k][16][t] = 0;
		} else if (k != 0) {
			AF_elem[m][k][16][t] = -(-(a_1 * cell_values[m][k - 1][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_1 *
							(((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 -
									(a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_2 *
							(((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 -
									(a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_3 *
							(((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 -
									(a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_4 *
							(((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				AF_elem[m][k][16][t] = -(-(-a_1 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_1 *
							(((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 -
									(a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_2 *
							(((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 -
									(a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_3 *
							(((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 -
									(a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_4 *
							(((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
			} else {
				AF_elem[m][k][16][t] = -(-(a_1 * cell_values[i_TOT - 2 - m][k][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_1 *
							(((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 -
									(a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_2 *
							(((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
							  (nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
							 ((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 -
									(a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_3 *
							(((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
							  (nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 -
									(a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t]) * dx_4 *
							(((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
							  (nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
							 ((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
							  (nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
			}
		}
	} else {
		AF_elem[m][k][16][t] = -( -dx_1 * ((aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t]) *
								((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
								(nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) +
								(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t]) *
								((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
								(nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])) * 0.5 -
						dx_2 * ((aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t]) *
								((nodal_values[m+1][k][0][t] - X_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t] -
								(nodal_values[m+1][k][2][t] - Z_POLE[t])*nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) +
								(aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t]) *
								((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
								(nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])) * 0.5 -
						dx_3 * ((aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t]) *
								((nodal_values[m+1][k+1][0][t] - X_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t] -
								(nodal_values[m+1][k+1][2][t] - Z_POLE[t])*nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) +
								(aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t]) *
								((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
								(nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])) * 0.5 -
						dx_4 * ((aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t]) *
								((nodal_values[m][k+1][0][t] - X_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t] -
								(nodal_values[m][k+1][2][t] - Z_POLE[t])*nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) +
								(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t]) *
								((nodal_values[m][k][0][t] - X_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t] -
								(nodal_values[m][k][2][t] - Z_POLE[t])*nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t])) * 0.5);
	}
	//---Second term
	if (flux_scheme) {
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			AF_elem[m][k][19][t] = 0;
		} else if (k != 0) {
			AF_elem[m][k][19][t]= (
				(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * cell_values[m][k - 1][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) *
											(-nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t] - nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) * 0.5 +
				(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) *
											(-nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t] - nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) * 0.5 +
				(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) *
											(-nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t] - nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) * 0.5 +
				(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t])*
											(-nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t] - nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) * 0.5);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				AF_elem[m][k][19][t]= (
					(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 + a_1) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(-nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t] - nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) * 0.5 +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(-nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t] - nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) * 0.5 +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(-nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t] - nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) * 0.5 +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t])*
												(-nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t] - nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) * 0.5);
			} else {
				AF_elem[m][k][19][t]= (
					(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * cell_values[i_TOT - 2 - m][k][numVars+compVarsN+12][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(-nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t] - nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t]) * 0.5 +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * cell_values[m + 1][k][numVars+compVarsN+12][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(-nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t] - nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]) * 0.5 +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * cell_values[m][k + 1][numVars+compVarsN+12][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+12][t]) *
												(-nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t] - nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]) * 0.5 +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * cell_values[m - 1][k][numVars+compVarsN+12][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+12][t])*
												(-nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t] - nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]) * 0.5);
			}
		}
	} else {
		AF_elem[m][k][19][t] = (
			(rx_1 * dz_1 - rz_1 * dx_1) * (-nodal_values[m][k][rho_id-1][t]*(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vx_id-1][t] - 
											nodal_values[m+1][k][rho_id-1][t]*(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vx_id-1][t]) * 0.5 +
			(rx_2 * dz_2 - rz_2 * dx_2) * (-nodal_values[m+1][k][rho_id-1][t]*(aux_values_N[m+1][k][22][t]-aux_values_N[m+1][k][23][t])*nodal_values[m+1][k][vx_id-1][t] - 
											nodal_values[m+1][k+1][rho_id-1][t]*(aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t])*nodal_values[m+1][k+1][vx_id-1][t]) * 0.5 +
			(rx_3 * dz_3 - rz_3 * dx_3) * (-nodal_values[m+1][k+1][rho_id-1][t]*(aux_values_N[m+1][k+1][22][t]-aux_values_N[m+1][k+1][23][t])*nodal_values[m+1][k+1][vx_id-1][t] -
											nodal_values[m][k+1][rho_id-1][t]*(aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t])*nodal_values[m][k+1][vx_id-1][t]) * 0.5 +
			(rx_4 * dz_4 - rz_4 * dx_4) * (-nodal_values[m][k+1][rho_id-1][t]*(aux_values_N[m][k+1][22][t]-aux_values_N[m][k+1][23][t])*nodal_values[m][k+1][vx_id-1][t] -
											nodal_values[m][k][rho_id-1][t]*(aux_values_N[m][k][22][t]-aux_values_N[m][k][23][t])*nodal_values[m][k][vx_id-1][t]) * 0.5);
	}
}