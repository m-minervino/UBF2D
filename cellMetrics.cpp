#include <math.h>
#include <vector>
#include "/opt/tecplot/360ex_2025r2/include/TECIO.h"

using namespace std;

void cell_metrics(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values, bool NS, int U_FORM, bool MOV_GRD, int t, int m, int k,
				  int numVars, int compVarsN, vector < bool > &near_ID, int i_TOT, int64_t* kMax, int grad_scheme, bool flux_scheme, bool ParForm, vector < double > &X_POLE, vector < double > &Z_POLE,
				  double& rx_1, double& rz_1, double& rx_2, double& rz_2, double& rx_3, double& rz_3, double& rx_4, double& rz_4,
				  double& dx_1, double& dz_1, double& dx_2, double& dz_2, double& dx_3, double& dz_3, double& dx_4, double& dz_4,
				  double& a_1, double& a_2, double& a_3, double& a_4, vector < vector < double > > &WLS) {
	double D[4][2], w[4], G[2][2], invG[2][2], detG=0;
	if (grad_scheme==2) {
		for (int i=0; i<4; i++) {	//Initialization loop
			w[i]=0;
			for (int j=0; j<2; j++) {
				D[i][j]=0;
				if (i<2) {
					G[i][j]=0;
					invG[i][j]=0;
				}
			}
		}
	}
	rx_1=((nodal_values[m+1][k][0][t]+nodal_values[m][k][0][t])/2)-X_POLE[t];
	rz_1=((nodal_values[m+1][k][2][t]+nodal_values[m][k][2][t])/2)-Z_POLE[t];
	rx_2=((nodal_values[m+1][k+1][0][t]+nodal_values[m+1][k][0][t])/2)-X_POLE[t];
	rz_2=((nodal_values[m+1][k+1][2][t]+nodal_values[m+1][k][2][t])/2)-Z_POLE[t];
	rx_3=((nodal_values[m][k+1][0][t]+nodal_values[m+1][k+1][0][t])/2)-X_POLE[t];
	rz_3=((nodal_values[m][k+1][2][t]+nodal_values[m+1][k+1][2][t])/2)-Z_POLE[t];
	rx_4=((nodal_values[m][k][0][t]+nodal_values[m][k+1][0][t])/2)-X_POLE[t];
	rz_4=((nodal_values[m][k][2][t]+nodal_values[m][k+1][2][t])/2)-Z_POLE[t];
	dx_1=nodal_values[m+1][k][0][t]-nodal_values[m][k][0][t];
	dz_1=nodal_values[m+1][k][2][t]-nodal_values[m][k][2][t];
	dx_2=nodal_values[m+1][k+1][0][t]-nodal_values[m+1][k][0][t];
	dz_2=nodal_values[m+1][k+1][2][t]-nodal_values[m+1][k][2][t];
	dx_3=nodal_values[m][k+1][0][t]-nodal_values[m+1][k+1][0][t];
	dz_3=nodal_values[m][k+1][2][t]-nodal_values[m+1][k+1][2][t];
	dx_4=nodal_values[m][k][0][t]-nodal_values[m][k+1][0][t];
	dz_4=nodal_values[m][k][2][t]-nodal_values[m][k+1][2][t];
	if (flux_scheme) {
		if (k != 0) {
			a_1 = sqrt(pow(rx_1 - cell_values[m][k][numVars + compVarsN + 15][t], 2) + pow(rz_1 - cell_values[m][k][numVars + compVarsN + 16][t], 2)) /
					sqrt(pow(cell_values[m][k - 1][numVars + compVarsN + 15][t] - cell_values[m][k][numVars + compVarsN + 15][t], 2) +
						pow(cell_values[m][k - 1][numVars + compVarsN + 16][t] - cell_values[m][k][numVars + compVarsN + 16][t], 2));
			//a_1 = sqrt(pow(rx_1 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_1 - cell_values[m][k][numVars + compVarsN + 16], 2)) /
			//		(sqrt(pow(rx_1 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_1 - cell_values[m][k][numVars + compVarsN + 16], 2))+
			//		sqrt(pow(rx_1 - cell_values[m][k-1][numVars + compVarsN + 15], 2) + pow(rz_1 - cell_values[m][k-1][numVars + compVarsN + 16], 2)));
		}
		else {
			if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
				if ( (!ParForm) && (NS && (!MOV_GRD || (U_FORM==1))) ) {
					a_1 = 0.5;
				} else {
					a_1 = sqrt(pow(rx_1 - cell_values[m][k][numVars + compVarsN + 15][t], 2) + pow(rz_1 - cell_values[m][k][numVars + compVarsN + 16][t], 2)) /
						  sqrt(pow(cell_values[m][k + 1][numVars + compVarsN + 15][t] - cell_values[m][k][numVars + compVarsN + 15][t], 2) +
							   pow(cell_values[m][k + 1][numVars + compVarsN + 16][t] - cell_values[m][k][numVars + compVarsN + 16][t], 2));
				}
			}
			else {
				a_1 = sqrt(pow(rx_1 - cell_values[m][k][numVars + compVarsN + 15][t], 2) + pow(rz_1 - cell_values[m][k][numVars + compVarsN + 16][t], 2)) /
						sqrt(pow(cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 15][t] - cell_values[m][k][numVars + compVarsN + 15][t], 2) +
							pow(cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 16][t] - cell_values[m][k][numVars + compVarsN + 16][t], 2));
			}
		}
		if (m != i_TOT - 2) {
			a_2 = sqrt(pow(rx_2 - cell_values[m][k][numVars + compVarsN + 15][t], 2) + pow(rz_2 - cell_values[m][k][numVars + compVarsN + 16][t], 2)) /
					sqrt(pow(cell_values[m + 1][k][numVars + compVarsN + 15][t] - cell_values[m][k][numVars + compVarsN + 15][t], 2) +
						pow(cell_values[m + 1][k][numVars + compVarsN + 16][t] - cell_values[m][k][numVars + compVarsN + 16][t], 2));
			//a_2 = sqrt(pow(rx_2 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_2 - cell_values[m][k][numVars + compVarsN + 16], 2)) /
			//		(sqrt(pow(rx_2 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_2 - cell_values[m][k][numVars + compVarsN + 16], 2))+
			//		sqrt(pow(rx_2 - cell_values[m+1][k][numVars + compVarsN + 15], 2) + pow(rz_2 - cell_values[m+1][k][numVars + compVarsN + 16], 2)));
		}
		else {
			a_2 = 0;
		}
		if (k != kMax[0] - 2) {
			a_3 = sqrt(pow(rx_3 - cell_values[m][k][numVars + compVarsN + 15][t], 2) + pow(rz_3 - cell_values[m][k][numVars + compVarsN + 16][t], 2)) /
					sqrt(pow(cell_values[m][k + 1][numVars + compVarsN + 15][t] - cell_values[m][k][numVars + compVarsN + 15][t], 2) +
						pow(cell_values[m][k + 1][numVars + compVarsN + 16][t] - cell_values[m][k][numVars + compVarsN + 16][t], 2));
			//a_3 = sqrt(pow(rx_3 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_3 - cell_values[m][k][numVars + compVarsN + 16], 2)) /
			//		(sqrt(pow(rx_3 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_3 - cell_values[m][k][numVars + compVarsN + 16], 2))+
			//		sqrt(pow(rx_3 - cell_values[m][k+1][numVars + compVarsN + 15], 2) + pow(rz_3 - cell_values[m][k+1][numVars + compVarsN + 16], 2)));
		}
		else {
			a_3 = 0;
		}
		if (m != 0) {
			a_4 = sqrt(pow(rx_4 - cell_values[m][k][numVars + compVarsN + 15][t], 2) + pow(rz_4 - cell_values[m][k][numVars + compVarsN + 16][t], 2)) /
					sqrt(pow(cell_values[m - 1][k][numVars + compVarsN + 15][t] - cell_values[m][k][numVars + compVarsN + 15][t], 2) +
						pow(cell_values[m - 1][k][numVars + compVarsN + 16][t] - cell_values[m][k][numVars + compVarsN + 16][t], 2));
			//a_4 = sqrt(pow(rx_4 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_4 - cell_values[m][k][numVars + compVarsN + 16], 2)) /
			//		(sqrt(pow(rx_4 - cell_values[m][k][numVars + compVarsN + 15], 2) + pow(rz_4 - cell_values[m][k][numVars + compVarsN + 16], 2))+
			//		sqrt(pow(rx_4 - cell_values[m-1][k][numVars + compVarsN + 15], 2) + pow(rz_4 - cell_values[m-1][k][numVars + compVarsN + 16], 2)));
		}
		else {
			a_4 = 0;
		}
	}
	if (grad_scheme==2) {
		//Assemble distance and weights matrices
			if (k != 0) {
				D[0][0]=cell_values[m][k - 1][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t];
				D[0][1]=cell_values[m][k - 1][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t];
			} else {
				if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
					// D[0][0]=-(cell_values[m][k + 1][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t]);
					// D[0][1]=-(cell_values[m][k + 1][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t]);
					//Do nothing. Retain zero initialization to use a three (instead of four) cells stencil.
				} else {
					D[0][0]=cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t];
					D[0][1]=cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t];
				}
			}
			if (m != i_TOT - 2) {
				D[1][0]=cell_values[m + 1][k][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t];
				D[1][1]=cell_values[m + 1][k][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t];
			} else {
				D[1][0]=-(cell_values[m - 1][k][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t]);
				D[1][1]=-(cell_values[m - 1][k][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t]);
			}
			if (k != kMax[0] - 2) {
				D[2][0]=cell_values[m][k + 1][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t];
				D[2][1]=cell_values[m][k + 1][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t];
			} else {
				D[2][0]=-(cell_values[m][k - 1][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t]);
				D[2][1]=-(cell_values[m][k - 1][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t]);
			}
			if (m != 0) {
				D[3][0]=cell_values[m - 1][k][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t];
				D[3][1]=cell_values[m - 1][k][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t];
			} else {
				D[3][0]=-(cell_values[m + 1][k][numVars + compVarsN + 15][t]-cell_values[m][k][numVars + compVarsN + 15][t]);
				D[3][1]=-(cell_values[m + 1][k][numVars + compVarsN + 16][t]-cell_values[m][k][numVars + compVarsN + 16][t]);
			}
			for (int i=0; i<4; i++) {
				if ((D[i][0]==0)&&(D[i][1]==0)) {
					w[i]=0;
				} else {
					w[i]=1/sqrt(pow(D[i][0], 2) + pow(D[i][1], 2));
				}
			}
		//Assemble G matrix
		G[0][0]=w[0]*w[0]*D[0][0]*D[0][0] + w[1]*w[1]*D[1][0]*D[1][0] + w[2]*w[2]*D[2][0]*D[2][0] + w[3]*w[3]*D[3][0]*D[3][0];
		G[0][1]=w[0]*w[0]*D[0][0]*D[0][1] + w[1]*w[1]*D[1][0]*D[1][1] + w[2]*w[2]*D[2][0]*D[2][1] + w[3]*w[3]*D[3][0]*D[3][1];
		G[1][0]=w[0]*w[0]*D[0][1]*D[0][0] + w[1]*w[1]*D[1][1]*D[1][0] + w[2]*w[2]*D[2][1]*D[2][0] + w[3]*w[3]*D[3][1]*D[3][0];
		G[1][1]=w[0]*w[0]*D[0][1]*D[0][1] + w[1]*w[1]*D[1][1]*D[1][1] + w[2]*w[2]*D[2][1]*D[2][1] + w[3]*w[3]*D[3][1]*D[3][1];
		//G matrix inverse
		detG=G[0][0]*G[1][1]-G[0][1]*G[1][0];
		invG[0][0]= G[1][1]/detG;
		invG[1][1]= G[0][0]/detG;
		invG[0][1]=-G[0][1]/detG;
		invG[1][0]=-G[1][0]/detG;
		//Assemble Weighted Least Squares matrix
		WLS[0][0]=w[0]*w[0]*(invG[0][0]*D[0][0] + invG[0][1]*D[0][1]);
		WLS[0][1]=w[1]*w[1]*(invG[0][0]*D[1][0] + invG[0][1]*D[1][1]);
		WLS[0][2]=w[2]*w[2]*(invG[0][0]*D[2][0] + invG[0][1]*D[2][1]);
		WLS[0][3]=w[3]*w[3]*(invG[0][0]*D[3][0] + invG[0][1]*D[3][1]);
		WLS[1][0]=w[0]*w[0]*(invG[1][0]*D[0][0] + invG[1][1]*D[0][1]);
		WLS[1][1]=w[1]*w[1]*(invG[1][0]*D[1][0] + invG[1][1]*D[1][1]);
		WLS[1][2]=w[2]*w[2]*(invG[1][0]*D[2][0] + invG[1][1]*D[2][1]);
		WLS[1][3]=w[3]*w[3]*(invG[1][0]*D[3][0] + invG[1][1]*D[3][1]);
	}
}