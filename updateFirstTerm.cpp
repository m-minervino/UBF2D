#include <vector>
#include "/opt/tecplot/360ex_2025r2/include/TECIO.h"

using namespace std;

void updateFirstTerm(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values, int rho_id,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, int U_FORM, bool MOV_GRD, int mrho_form, bool flux_scheme, vector < bool > &near_ID,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem) {
	//******** WARNING *********
	//This function computes curl (mrho_form==0 and mrho_form==1) using the G-G nodes-based method only.
	//Other gradient schemes are to be included and implemented, as done in the rest of the code.
	//**************************
	if (mrho_form == 0) {
		//Curl of rho*grad(0.5*V^2) times cell volume
		double curl=0;
		if (flux_scheme) {
			if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
				// curl is set at zero as far-field cells will never be used. Therefore no action is needed as curl initialization is enough.
			} else if (k != 0) {
				curl = (0.5*(nodal_values[m][k][rho_id - 1][t]+nodal_values[m+1][k][rho_id - 1][t])*(
						(-(a_1 * cell_values[m][k-1][numVars + compVarsN + 11][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_1)+
						(-(a_1 * cell_values[m][k-1][numVars + compVarsN + 10][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_1))+
						0.5*(nodal_values[m+1][k][rho_id - 1][t]+nodal_values[m+1][k+1][rho_id - 1][t])*(
						(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 11][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_2)+
						(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 10][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_2))+
						0.5*(nodal_values[m+1][k+1][rho_id - 1][t]+nodal_values[m][k+1][rho_id - 1][t])*(
						(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 11][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_3)+
						(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 10][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_3))+
						0.5*(nodal_values[m][k+1][rho_id - 1][t]+nodal_values[m][k][rho_id - 1][t])*(
						(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 11][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_4)+
						(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 10][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_4)));
			} else {
				if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
					if (( NS && (!MOV_GRD || (U_FORM==1)))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option
						curl = (0.5*(nodal_values[m][k][rho_id - 1][t]+nodal_values[m+1][k][rho_id - 1][t])*(
								(-(a_1 * cell_values[m][k][numVars + compVarsN + 11][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_1)+
								(-(a_1 * cell_values[m][k][numVars + compVarsN + 10][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_1))+
								0.5*(nodal_values[m+1][k][rho_id - 1][t]+nodal_values[m+1][k+1][rho_id - 1][t])*(
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 11][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_2)+
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 10][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_2))+
								0.5*(nodal_values[m+1][k+1][rho_id - 1][t]+nodal_values[m][k+1][rho_id - 1][t])*(
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 11][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_3)+
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 10][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_3))+
								0.5*(nodal_values[m][k+1][rho_id - 1][t]+nodal_values[m][k][rho_id - 1][t])*(
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 11][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_4)+
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 10][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_4)));
					} else {	//In this case the cell_metrics function has set a_1 to allow extrapolation based on the first two cells.
						curl = (0.5*(nodal_values[m][k][rho_id - 1][t]+nodal_values[m+1][k][rho_id - 1][t])*(
								(-(-a_1 * cell_values[m][k+1][numVars + compVarsN + 11][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_1)+
								(-(-a_1 * cell_values[m][k+1][numVars + compVarsN + 10][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_1))+
								0.5*(nodal_values[m+1][k][rho_id - 1][t]+nodal_values[m+1][k+1][rho_id - 1][t])*(
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 11][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_2)+
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 10][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_2))+
								0.5*(nodal_values[m+1][k+1][rho_id - 1][t]+nodal_values[m][k+1][rho_id - 1][t])*(
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 11][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_3)+
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 10][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_3))+
								0.5*(nodal_values[m][k+1][rho_id - 1][t]+nodal_values[m][k][rho_id - 1][t])*(
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 11][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_4)+
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 10][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_4)));
					}
				} else {
					curl = (0.5*(nodal_values[m][k][rho_id - 1][t]+nodal_values[m+1][k][rho_id - 1][t])*(
							(-(a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 11][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_1)+
							(-(a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 10][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_1))+
							0.5*(nodal_values[m+1][k][rho_id - 1][t]+nodal_values[m+1][k+1][rho_id - 1][t])*(
							(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 11][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_2)+
							(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 10][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_2))+
							0.5*(nodal_values[m+1][k+1][rho_id - 1][t]+nodal_values[m][k+1][rho_id - 1][t])*(
							(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 11][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_3)+
							(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 10][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_3))+
							0.5*(nodal_values[m][k+1][rho_id - 1][t]+nodal_values[m][k][rho_id - 1][t])*(
							(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 11][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 11][t])*dz_4)+
							(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 10][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 10][t])*dx_4)));
				}
			}
		} else {
			// To be implemented
		}
		AF_elem[m][k][4][t] = cell_values[m][k][numVars + compVarsN + 16][t] * curl;
		AF_elem[m][k][5][t] = cell_values[m][k][numVars + compVarsN + 15][t] * (-curl);
	} else if (mrho_form == 1) {
		//Curl of 0.5*V^2*grad(rho)
		double curl=0;
		if (flux_scheme) {
			if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
				// curl is set at zero as far-field cells will never be used. Therefore no action is needed as curl initialization is enough.
			} else if (k != 0) {
				curl = (0.5*(nodal_values[m][k][numVars + 1][t]+nodal_values[m+1][k][numVars + 1][t])*(
						(-(a_1 * cell_values[m][k-1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_1)+
						(-(a_1 * cell_values[m][k-1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_1))+
						0.5*(nodal_values[m+1][k][numVars + 1][t]+nodal_values[m+1][k+1][numVars + 1][t])*(
						(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_2)+
						(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_2))+
						0.5*(nodal_values[m+1][k+1][numVars + 1][t]+nodal_values[m][k+1][numVars + 1][t])*(
						(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_3)+
						(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_3))+
						0.5*(nodal_values[m][k+1][numVars + 1][t]+nodal_values[m][k][numVars + 1][t])*(
						(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_4)+
						(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_4)));
			} else {
				if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
					if (( NS && (!MOV_GRD || (U_FORM==1)))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option
						curl = (0.5*(nodal_values[m][k][numVars + 1][t]+nodal_values[m+1][k][numVars + 1][t])*(
								(-(a_1 * cell_values[m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_1)+
								(-(a_1 * cell_values[m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_1))+
								0.5*(nodal_values[m+1][k][numVars + 1][t]+nodal_values[m+1][k+1][numVars + 1][t])*(
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_2)+
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_2))+
								0.5*(nodal_values[m+1][k+1][numVars + 1][t]+nodal_values[m][k+1][numVars + 1][t])*(
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_3)+
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_3))+
								0.5*(nodal_values[m][k+1][numVars + 1][t]+nodal_values[m][k][numVars + 1][t])*(
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_4)+
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_4)));
					} else {	//In this case the cell_metrics function has set a_1 to allow extrapolation based on the first two cells.
						curl = (0.5*(nodal_values[m][k][numVars + 1][t]+nodal_values[m+1][k][numVars + 1][t])*(
								(-(-a_1 * cell_values[m][k+1][numVars + compVarsN + 9][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_1)+
								(-(-a_1 * cell_values[m][k+1][numVars + compVarsN + 8][t] + (1 + a_1) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_1))+
								0.5*(nodal_values[m+1][k][numVars + 1][t]+nodal_values[m+1][k+1][numVars + 1][t])*(
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_2)+
								(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_2))+
								0.5*(nodal_values[m+1][k+1][numVars + 1][t]+nodal_values[m][k+1][numVars + 1][t])*(
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_3)+
								(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_3))+
								0.5*(nodal_values[m][k+1][numVars + 1][t]+nodal_values[m][k][numVars + 1][t])*(
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_4)+
								(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_4)));
					}
				} else {
					curl = (0.5*(nodal_values[m][k][numVars + 1][t]+nodal_values[m+1][k][numVars + 1][t])*(
							(-(a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_1)+
							(-(a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_1))+
							0.5*(nodal_values[m+1][k][numVars + 1][t]+nodal_values[m+1][k+1][numVars + 1][t])*(
							(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_2)+
							(-(a_2 * cell_values[m+1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_2))+
							0.5*(nodal_values[m+1][k+1][numVars + 1][t]+nodal_values[m][k+1][numVars + 1][t])*(
							(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_3)+
							(-(a_3 * cell_values[m][k+1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_3))+
							0.5*(nodal_values[m][k+1][numVars + 1][t]+nodal_values[m][k][numVars + 1][t])*(
							(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t])*dz_4)+
							(-(a_4 * cell_values[m-1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t])*dx_4)));
				}
			}
		} else {
			// To be implemented
		}
		AF_elem[m][k][4][t] = cell_values[m][k][numVars + compVarsN + 16][t] * (-curl);
		AF_elem[m][k][5][t] = cell_values[m][k][numVars + compVarsN + 15][t] * curl;
	} else if (mrho_form==2) {
		AF_elem[m][k][4][t] = -cell_values[m][k][numVars + compVarsN + 16][t] * (cell_values[m][k][numVars + compVarsN + 8][t] * cell_values[m][k][numVars + compVarsN + 11][t] -
			cell_values[m][k][numVars + compVarsN + 9][t] * cell_values[m][k][numVars + compVarsN + 10][t]) * cell_values[m][k][numVars + compVarsN][t];
		AF_elem[m][k][5][t] = -cell_values[m][k][numVars + compVarsN + 15][t] * (-cell_values[m][k][numVars + compVarsN + 8][t] * cell_values[m][k][numVars + compVarsN + 11][t] +
			cell_values[m][k][numVars + compVarsN + 9][t] * cell_values[m][k][numVars + compVarsN + 10][t]) * cell_values[m][k][numVars + compVarsN][t];
	}
}