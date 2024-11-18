#include <vector>

using namespace std;

void updateFourthTerm(vector < vector < vector < vector < double > > > > &nodal_values,vector < vector < vector < vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < double > > > > &cell_values,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, bool flux_scheme, vector < bool > &near_ID,
					double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < bool > > > > &flags,
					vector < vector < vector < vector < double > > > > &AF_elem) {
	if (k != 0) {
		if (flags[m][k - 1][0][t] == false) {
			if (flux_scheme) {
				AF_elem[m][k][10][t] = AF_elem[m][k][10][t] +
					rz_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
						dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
						dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
			}
			else {
				AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + 
					rz_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
						dz_1 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m][k][18][t]) / 2) +
						dx_1 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m][k][17][t]) / 2));
			}
		}
	}
	else {
		if ((flags[i_TOT - 2 - m][k][0][t] == false) && (!near_ID[m]) && (!near_ID[m + 1])) {
			if (flux_scheme) {
				AF_elem[m][k][10][t] = AF_elem[m][k][10][t] +
					rz_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
						dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
						dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
			}
			else {
				AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + 
					rz_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
						dz_1 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m][k][18][t]) / 2) +
						dx_1 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m][k][17][t]) / 2));
			}
		} else if ((near_ID[m]) || (near_ID[m + 1])) {
			//Not effective in case of NS and fixed grid computation as V will be zero on the body surface
			if (flux_scheme) {
				AF_elem[m][k][10][t] = AF_elem[m][k][10][t] +
					rz_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
						dz_1 * (-a_1 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
						dx_1 * (-a_1 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
			}
			else {
				AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + 
					rz_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
						dz_1 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m][k][18][t]) / 2) +
						dx_1 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m][k][17][t]) / 2));
			}
		}
	}
	if (flags[m+1][k][0][t]==false) {
		if (flux_scheme) {
			AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + rz_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
		} else {
			AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + rz_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
													dx_2 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
		}
	}
	if (flags[m][k+1][0][t]==false) {
		if (flux_scheme) {
			AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
		} else {
			AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + rz_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
													dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
		}
	}
	if (flags[m-1][k][0][t]==false) {
		if (flux_scheme) {
			AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
		} else {
			AF_elem[m][k][10][t] = AF_elem[m][k][10][t] + rz_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
													dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
		}
	}
	if (k != 0) {
		if (flags[m][k - 1][0][t] == false) {
			if (flux_scheme) {
				AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - 
					rx_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
					dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
					dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
			}
			else {
				AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - 
					rx_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
					dz_1 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m][k][18][t]) / 2) +
					dx_1 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m][k][17][t]) / 2));
			}
		}
	}
	else {
		if ((flags[i_TOT - 2 - m][k][0][t] == false) && (!near_ID[m]) && (!near_ID[m + 1])) {
			if (flux_scheme) {
				AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - 
					rx_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
					dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
					dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
			}
			else {
				AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - 
					rx_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
					dz_1 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m][k][18][t]) / 2) +
					dx_1 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m][k][17][t]) / 2));
			}
		} else if ((near_ID[m]) || (near_ID[m + 1])) {
			//Not effective in case of NS and fixed grid computation as V will be zero on the body surface
			if (flux_scheme) {
				AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - 
					rx_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
					dz_1 * (-a_1 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
					dx_1 * (-a_1 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
			}
			else {
				AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - 
					rx_1 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) * (
					dz_1 * ((aux_values_N[m + 1][k][18][t] + aux_values_N[m][k][18][t]) / 2) +
					dx_1 * ((aux_values_N[m + 1][k][17][t] + aux_values_N[m][k][17][t]) / 2));
			}
		}
	}
	if (flags[m+1][k][0][t]==false) {
		if (flux_scheme) {
			AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - rx_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
													(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
		} else {
			AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - rx_2 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
													(dz_2 * ((aux_values_N[m + 1][k + 1][18][t] + aux_values_N[m + 1][k][18][t]) / 2) +
													dx_2 * ((aux_values_N[m + 1][k + 1][17][t] + aux_values_N[m + 1][k][17][t]) / 2));
		}
	}
	if (flags[m][k+1][0][t]==false) {
		if (flux_scheme) {
			AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
		} else {
			AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - rx_3 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
													(dz_3 * ((aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 2) +
													dx_3 * ((aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 2));
		}
	}
	if (flags[m-1][k][0][t]==false) {
		if (flux_scheme) {
			AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
		} else {
			AF_elem[m][k][11][t] = AF_elem[m][k][11][t] - rx_4 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
													(dz_4 * ((aux_values_N[m][k][18][t] + aux_values_N[m][k + 1][18][t]) / 2) +
													dx_4 * ((aux_values_N[m][k][17][t] + aux_values_N[m][k + 1][17][t]) / 2));
		}
	}
}