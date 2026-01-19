#include <vector>
#include <math.h>
#include "/opt/tecplot/360ex_2025r2/include/TECIO.h"

using namespace std;

void updateF_grad_rho(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values,
					vector < vector < vector < vector < double > > > > &aux_values, vector < vector < vector < vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < bool > > > > &flags,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, int U_FORM, bool MOV_GRD, bool flux_scheme, bool smart_wds,
					vector < bool > &near_ID, double CdP_tsh, double rho_inf, double V_inf,
					double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem) {
	//******** WARNING *********
	//This function computes gradients (and divergence) using the G-G nodes-based method only.
	//Other gradient schemes are to be included and implemented, as done in the rest of the code.
	//**************************
	//---X-component
	double CdeP_vrt_gradrho=0;
	//------Density gradient term
	CdeP_vrt_gradrho = CdeP_vrt_gradrho + 0.5*pow(V_inf,2)*cell_values[m][k][numVars+compVarsN+8][t]*cell_values[m][k][numVars+compVarsN][t];
	//------Other terms
	if (flux_scheme) {
		//------First term
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			CdeP_vrt_gradrho = CdeP_vrt_gradrho + 0;
		} else if (k != 0) {
			CdeP_vrt_gradrho = CdeP_vrt_gradrho + ((a_1 * aux_values[m][k - 1][9][t] + (1 - a_1) * aux_values[m][k][9][t]) * dz_1 +
												   (a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dz_2 +
												   (a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dz_3 +
												   (a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dz_4);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || (U_FORM==1)))) {
					CdeP_vrt_gradrho = CdeP_vrt_gradrho + ((-a_1 * aux_values[m][k][9][t] + (1 - a_1) * aux_values[m][k][9][t]) * dz_1 +
														   (a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dz_2 +
														   (a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dz_3 +
														   (a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dz_4);
				} else {
					CdeP_vrt_gradrho = CdeP_vrt_gradrho + ((-a_1 * aux_values[m][k + 1][9][t] + (1 + a_1) * aux_values[m][k][9][t]) * dz_1 +
														   (a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dz_2 +
														   (a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dz_3 +
														   (a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dz_4);
				}
			} else {
				CdeP_vrt_gradrho = CdeP_vrt_gradrho + ((a_1 * aux_values[i_TOT - 2 - m][k][9][t] + (1 - a_1) * aux_values[m][k][9][t]) * dz_1 +
													   (a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dz_2 +
													   (a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dz_3 +
													   (a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dz_4);
			}
		}
		//------Second term
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			CdeP_vrt_gradrho = CdeP_vrt_gradrho - 0;
		} else if (k != 0) {
			CdeP_vrt_gradrho = CdeP_vrt_gradrho - (
				(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[m][k - 1][10][t] + (1 - a_1) * aux_values[m][k][10][t]) +
				(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][10][t] + (1 - a_2) * aux_values[m][k][10][t]) +
				(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][10][t] + (1 - a_3) * aux_values[m][k][10][t]) +
				(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][10][t] + (1 - a_4) * aux_values[m][k][10][t]));
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || (U_FORM==1)))) {
					CdeP_vrt_gradrho = CdeP_vrt_gradrho - (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k][10][t] + (1 - a_1) * aux_values[m][k][10][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][10][t] + (1 - a_2) * aux_values[m][k][10][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][10][t] + (1 - a_3) * aux_values[m][k][10][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][10][t] + (1 - a_4) * aux_values[m][k][10][t]));
				} else {
					CdeP_vrt_gradrho = CdeP_vrt_gradrho - (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k + 1][10][t] + (1 + a_1) * aux_values[m][k][10][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][10][t] + (1 - a_2) * aux_values[m][k][10][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][10][t] + (1 - a_3) * aux_values[m][k][10][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][10][t] + (1 - a_4) * aux_values[m][k][10][t]));
				}
			} else {
				CdeP_vrt_gradrho = CdeP_vrt_gradrho - (
					(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[i_TOT - 2 - m][k][10][t] + (1 - a_1) * aux_values[m][k][10][t]) +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][10][t] + (1 - a_2) * aux_values[m][k][10][t]) +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][10][t] + (1 - a_3) * aux_values[m][k][10][t]) +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][10][t] + (1 - a_4) * aux_values[m][k][10][t]));
			}
		}
	} else {
		//------First term
		CdeP_vrt_gradrho = CdeP_vrt_gradrho + (((aux_values_N[m][k][14][t] + aux_values_N[m + 1][k][14][t]) / 2) * dz_1 +
											   ((aux_values_N[m + 1][k][14][t] + aux_values_N[m + 1][k + 1][14][t]) / 2) * dz_2 +
											   ((aux_values_N[m][k + 1][14][t] + aux_values_N[m + 1][k + 1][14][t]) / 2) * dz_3 +
											   ((aux_values_N[m][k][14][t] + aux_values_N[m][k + 1][14][t]) / 2) * dz_4);
		//------Second term
		CdeP_vrt_gradrho = CdeP_vrt_gradrho - (
			(rx_1 * dz_1 - rz_1 * dx_1) * ((aux_values_N[m][k][15][t] + aux_values_N[m + 1][k][15][t]) / 2) +
			(rx_2 * dz_2 - rz_2 * dx_2) * ((aux_values_N[m + 1][k][15][t] + aux_values_N[m + 1][k + 1][15][t]) / 2) +
			(rx_3 * dz_3 - rz_3 * dx_3) * ((aux_values_N[m][k + 1][15][t] + aux_values_N[m + 1][k + 1][15][t]) / 2) +
			(rx_4 * dz_4 - rz_4 * dx_4) * ((aux_values_N[m][k][15][t] + aux_values_N[m][k + 1][15][t]) / 2));
	}
	//X-component attribution
	if (smart_wds) {
		
	} else {
		if (flags[m][k][2][t]) {
			AF_elem[m][k][22][t] = CdeP_vrt_gradrho;
		} else if (flags[m][k][1][t]) {
			AF_elem[m][k][23][t] = CdeP_vrt_gradrho;
		} else {
			AF_elem[m][k][24][t] = CdeP_vrt_gradrho;
		}
	}
	//---Z-component
	//------Density gradient term
	AF_elem[m][k][25][t] = AF_elem[m][k][25][t] + 0.5*pow(V_inf,2)*cell_values[m][k][numVars+compVarsN+9][t]*cell_values[m][k][numVars+compVarsN][t];
	//------Other terms
	if (flux_scheme) {
		//------First term
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			AF_elem[m][k][25][t] = AF_elem[m][k][25][t] + 0;
		} else if (k != 0) {
			AF_elem[m][k][25][t] = AF_elem[m][k][25][t] = + (-(a_1 * aux_values[m][k - 1][9][t] + (1 - a_1) * aux_values[m][k][9][t]) * dx_1 -
														(a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dx_2 -
														(a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dx_3 -
														(a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dx_4);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || (U_FORM==1)))) {
					AF_elem[m][k][25][t] = AF_elem[m][k][25][t] = + (-(-a_1 * aux_values[m][k][9][t] + (1 - a_1) * aux_values[m][k][9][t]) * dx_1 -
																(a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dx_2 -
																(a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dx_3 -
																(a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dx_4);
				} else {
					AF_elem[m][k][25][t] = AF_elem[m][k][25][t] = + (-(-a_1 * aux_values[m][k + 1][9][t] + (1 + a_1) * aux_values[m][k][9][t]) * dx_1 -
																(a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dx_2 -
																(a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dx_3 -
																(a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dx_4);
				}
			} else {
				AF_elem[m][k][25][t] = AF_elem[m][k][25][t] = + (-(a_1 * aux_values[i_TOT - 2 - m][k][9][t] + (1 - a_1) * aux_values[m][k][9][t]) * dx_1 -
															(a_2 * aux_values[m + 1][k][9][t] + (1 - a_2) * aux_values[m][k][9][t]) * dx_2 -
															(a_3 * aux_values[m][k + 1][9][t] + (1 - a_3) * aux_values[m][k][9][t]) * dx_3 -
															(a_4 * aux_values[m - 1][k][9][t] + (1 - a_4) * aux_values[m][k][9][t]) * dx_4);
			}
		}
		//------Second term
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			AF_elem[m][k][25][t] = AF_elem[m][k][25][t] - 0;
		} else if (k != 0) {
			AF_elem[m][k][25][t] = AF_elem[m][k][25][t] - (
				(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[m][k - 1][11][t] + (1 - a_1) * aux_values[m][k][11][t]) +
				(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][11][t] + (1 - a_2) * aux_values[m][k][11][t]) +
				(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][11][t] + (1 - a_3) * aux_values[m][k][11][t]) +
				(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][11][t] + (1 - a_4) * aux_values[m][k][11][t]));
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || (U_FORM==1)))) {
					AF_elem[m][k][25][t] = AF_elem[m][k][25][t] - (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k][11][t] + (1 - a_1) * aux_values[m][k][11][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][11][t] + (1 - a_2) * aux_values[m][k][11][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][11][t] + (1 - a_3) * aux_values[m][k][11][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][11][t] + (1 - a_4) * aux_values[m][k][11][t]));
				} else {
					AF_elem[m][k][25][t] = AF_elem[m][k][25][t] - (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k + 1][11][t] + (1 + a_1) * aux_values[m][k][11][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][11][t] + (1 - a_2) * aux_values[m][k][11][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][11][t] + (1 - a_3) * aux_values[m][k][11][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][11][t] + (1 - a_4) * aux_values[m][k][11][t]));
				}
			} else {
				AF_elem[m][k][25][t] = AF_elem[m][k][25][t] - (
					(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[i_TOT - 2 - m][k][11][t] + (1 - a_1) * aux_values[m][k][11][t]) +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][11][t] + (1 - a_2) * aux_values[m][k][11][t]) +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][11][t] + (1 - a_3) * aux_values[m][k][11][t]) +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][11][t] + (1 - a_4) * aux_values[m][k][11][t]));
			}
		}
	} else {
		//------First term
		AF_elem[m][k][25][t] = AF_elem[m][k][25][t] + (-((aux_values_N[m][k][14][t] + aux_values_N[m + 1][k][14][t]) / 2) * dx_1 -
												  ((aux_values_N[m + 1][k][14][t] + aux_values_N[m + 1][k + 1][14][t]) / 2) * dx_2 -
												  ((aux_values_N[m][k + 1][14][t] + aux_values_N[m + 1][k + 1][14][t]) / 2) * dx_3 -
												  ((aux_values_N[m][k][14][t] + aux_values_N[m][k + 1][14][t]) / 2) * dx_4);
		//------Second term
		AF_elem[m][k][25][t] = AF_elem[m][k][25][t] - (
			(rx_1 * dz_1 - rz_1 * dx_1) * ((aux_values_N[m][k][16][t] + aux_values_N[m + 1][k][16][t]) / 2) +
			(rx_2 * dz_2 - rz_2 * dx_2) * ((aux_values_N[m + 1][k][16][t] + aux_values_N[m + 1][k + 1][16][t]) / 2) +
			(rx_3 * dz_3 - rz_3 * dx_3) * ((aux_values_N[m][k + 1][16][t] + aux_values_N[m + 1][k + 1][16][t]) / 2) +
			(rx_4 * dz_4 - rz_4 * dx_4) * ((aux_values_N[m][k][16][t] + aux_values_N[m][k + 1][16][t]) / 2));
	}
}