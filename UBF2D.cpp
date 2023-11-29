#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <chrono>

//Required by TecIO
#include "/usr/local/apps/tecplot2022r1/360ex_2022r1/include/TECIO.h"
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

//Linear Algebra (MTL4 library)
#include </u1/mima619/MTL4_GITHUB/trunk/boost/numeric/mtl/mtl.hpp>
#include </u1/mima619/MTL4_GITHUB/trunk/boost/numeric/itl/itl.hpp>

using namespace std;
using namespace mtl;
using namespace itl;

double S_tria(double x1, double z1, double x2, double z2, double x3, double z3) {
	return (-(x3-x1)*(z2-z1)+(z3-z1)*(x2-x1))/2;
}

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
				if ( (!ParForm) && (NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2)))) ) {
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

void updateCompCorr(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values, int rho_id,
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
					if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option
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
		AF_elem[m][k][2][t] = cell_values[m][k][numVars + compVarsN + 16][t] * curl;
		AF_elem[m][k][3][t] = cell_values[m][k][numVars + compVarsN + 15][t] * (-curl);
	} else if (mrho_form == 1) {
		//Curl of 0.5*V^2*grad(rho) times cell volume
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
					if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option
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
		AF_elem[m][k][2][t] = cell_values[m][k][numVars + compVarsN + 16][t] * (-curl);
		AF_elem[m][k][3][t] = cell_values[m][k][numVars + compVarsN + 15][t] * curl;
	} else if (mrho_form == 2) {
		AF_elem[m][k][2][t] = -cell_values[m][k][numVars + compVarsN + 16][t] * (cell_values[m][k][numVars + compVarsN + 8][t] * cell_values[m][k][numVars + compVarsN + 11][t] -
			cell_values[m][k][numVars + compVarsN + 9][t] * cell_values[m][k][numVars + compVarsN + 10][t]) * cell_values[m][k][numVars + compVarsN][t];
		AF_elem[m][k][3][t] = -cell_values[m][k][numVars + compVarsN + 15][t] * (-cell_values[m][k][numVars + compVarsN + 8][t] * cell_values[m][k][numVars + compVarsN + 11][t] +
			cell_values[m][k][numVars + compVarsN + 9][t] * cell_values[m][k][numVars + compVarsN + 10][t]) * cell_values[m][k][numVars + compVarsN][t];
	}
}

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
					if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option
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
					if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option
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

void updateParasite(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector <vector < double > > > > &cell_values,
					vector < vector < vector < vector < double > > > > &aux_values, vector < vector < vector <vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < bool > > > > &flags,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, int U_FORM, bool MOV_GRD, bool flux_scheme, bool smart_wds,
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
	//------First term
	if (flux_scheme) {	
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			CdeP_vrt_1 = 0;
		} else if (k != 0) {
			CdeP_vrt_1 = -((a_1 * aux_values[m][k - 1][0][t] + (1 - a_1) * aux_values[m][k][0][t]) * dz_1 +
									(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dz_2 +
									(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dz_3 +
									(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dz_4);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
					CdeP_vrt_1 = -((-a_1 * aux_values[m][k][0][t] + (1 - a_1) * aux_values[m][k][0][t]) * dz_1 +
									(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dz_2 +
									(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dz_3 +
									(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dz_4);
				} else {
					CdeP_vrt_1 = -((-a_1 * aux_values[m][k + 1][0][t] + (1 + a_1) * aux_values[m][k][0][t]) * dz_1 +
									(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dz_2 +
									(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dz_3 +
									(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dz_4);
				}
			} else {
				CdeP_vrt_1 = -((a_1 * aux_values[i_TOT - 2 - m][k][0][t] + (1 - a_1) * aux_values[m][k][0][t]) * dz_1 +
								(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dz_2 +
								(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dz_3 +
								(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dz_4);
			}
		}
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			CdeP_vrt_2 = 0;
		} else if (k != 0) {
			CdeP_vrt_2 = (
				(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[m][k - 1][1][t] + (1 - a_1) * aux_values[m][k][1][t]) +
				(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][1][t] + (1 - a_2) * aux_values[m][k][1][t]) +
				(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][1][t] + (1 - a_3) * aux_values[m][k][1][t]) +
				(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][1][t] + (1 - a_4) * aux_values[m][k][1][t]));
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
					CdeP_vrt_2 = (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k][1][t] + (1 - a_1) * aux_values[m][k][1][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][1][t] + (1 - a_2) * aux_values[m][k][1][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][1][t] + (1 - a_3) * aux_values[m][k][1][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][1][t] + (1 - a_4) * aux_values[m][k][1][t]));
				} else {
					CdeP_vrt_2 = (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k + 1][1][t] + (1 + a_1) * aux_values[m][k][1][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][1][t] + (1 - a_2) * aux_values[m][k][1][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][1][t] + (1 - a_3) * aux_values[m][k][1][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][1][t] + (1 - a_4) * aux_values[m][k][1][t]));
				}
			} else {
				CdeP_vrt_2 = (
					(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[i_TOT - 2 - m][k][1][t] + (1 - a_1) * aux_values[m][k][1][t]) +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][1][t] + (1 - a_2) * aux_values[m][k][1][t]) +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][1][t] + (1 - a_3) * aux_values[m][k][1][t]) +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][1][t] + (1 - a_4) * aux_values[m][k][1][t]));
			}
		}
		if (smart_wds) {
			if ((((CdeP_vrt_1+ CdeP_vrt_2) > CdP_tsh* 0.5 * rho_inf * V_inf * V_inf * chord) || ((CdeP_vrt_1 + CdeP_vrt_2) < -CdP_tsh * 0.5 * rho_inf * V_inf * V_inf * chord)) && flags[m][k][2][t]) {
				AF_elem[m][k][14][t] = CdeP_vrt_1;
			} else if (flags[m][k][1][t]) {
				AF_elem[m][k][15][t] = CdeP_vrt_1;
			} else {
				AF_elem[m][k][33][t] = CdeP_vrt_1;
				flags[m][k][2][t] = false;
			}
		} else {
			if (flags[m][k][2][t]) {
				AF_elem[m][k][14][t] = CdeP_vrt_1;
			} else if (flags[m][k][1][t]) {
				AF_elem[m][k][15][t] = CdeP_vrt_1;
			} else {
				AF_elem[m][k][33][t] = CdeP_vrt_1;
			}
		}
	} else {
		CdeP_vrt_1= -(((aux_values_N[m][k][0][t] + aux_values_N[m + 1][k][0][t]) / 2) * dz_1 +
			((aux_values_N[m + 1][k][0][t] + aux_values_N[m + 1][k + 1][0][t]) / 2) * dz_2 +
			((aux_values_N[m][k + 1][0][t] + aux_values_N[m + 1][k + 1][0][t]) / 2) * dz_3 +
			((aux_values_N[m][k][0][t] + aux_values_N[m][k + 1][0][t]) / 2) * dz_4);
		CdeP_vrt_2 = (
			(rx_1 * dz_1 - rz_1 * dx_1) * ((aux_values_N[m][k][1][t] + aux_values_N[m + 1][k][1][t]) / 2) +
			(rx_2 * dz_2 - rz_2 * dx_2) * ((aux_values_N[m + 1][k][1][t] + aux_values_N[m + 1][k + 1][1][t]) / 2) +
			(rx_3 * dz_3 - rz_3 * dx_3) * ((aux_values_N[m][k + 1][1][t] + aux_values_N[m + 1][k + 1][1][t]) / 2) +
			(rx_4 * dz_4 - rz_4 * dx_4) * ((aux_values_N[m][k][1][t] + aux_values_N[m][k + 1][1][t]) / 2));
		if (smart_wds) {
			if ((((CdeP_vrt_1+ CdeP_vrt_2) > CdP_tsh* 0.5 * rho_inf * V_inf * V_inf * chord) || ((CdeP_vrt_1 + CdeP_vrt_2) < -CdP_tsh * 0.5 * rho_inf * V_inf * V_inf * chord)) && flags[m][k][2][t]) {
				AF_elem[m][k][14][t] = CdeP_vrt_1;
			} else if (flags[m][k][1][t]) {
				AF_elem[m][k][15][t] = CdeP_vrt_1;
			} else {
				AF_elem[m][k][33][t] = CdeP_vrt_1;
				flags[m][k][2][t] = false;
			}
		} else {
			if (flags[m][k][2][t]) {
				AF_elem[m][k][14][t] = CdeP_vrt_1;
			} else if (flags[m][k][1][t]) {
				AF_elem[m][k][15][t] = CdeP_vrt_1;
			} else {
				AF_elem[m][k][33][t] = CdeP_vrt_1;
			}
		}
	}
	if (flux_scheme) {
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			AF_elem[m][k][16][t] = 0;
		} else if (k != 0) {
			AF_elem[m][k][16][t] = -(-(a_1 * aux_values[m][k - 1][0][t] + (1 - a_1) * aux_values[m][k][0][t]) * dx_1 -
									(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dx_2 -
									(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dx_3 -
									(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dx_4);
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
					AF_elem[m][k][16][t] = -(-(-a_1 * aux_values[m][k][0][t] + (1 - a_1) * aux_values[m][k][0][t]) * dx_1 -
											(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dx_2 -
											(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dx_3 -
											(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dx_4);
				} else {
					AF_elem[m][k][16][t] = -(-(-a_1 * aux_values[m][k + 1][0][t] + (1 + a_1) * aux_values[m][k][0][t]) * dx_1 -
											(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dx_2 -
											(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dx_3 -
											(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dx_4);
				}
			} else {
				AF_elem[m][k][16][t] = -(-(a_1 * aux_values[i_TOT - 2 - m][k][0][t] + (1 - a_1) * aux_values[m][k][0][t]) * dx_1 -
										(a_2 * aux_values[m + 1][k][0][t] + (1 - a_2) * aux_values[m][k][0][t]) * dx_2 -
										(a_3 * aux_values[m][k + 1][0][t] + (1 - a_3) * aux_values[m][k][0][t]) * dx_3 -
										(a_4 * aux_values[m - 1][k][0][t] + (1 - a_4) * aux_values[m][k][0][t]) * dx_4);
			}
		}
	} else {
		AF_elem[m][k][16][t] = -(-((aux_values_N[m][k][0][t] + aux_values_N[m + 1][k][0][t]) / 2) * dx_1 -
								((aux_values_N[m + 1][k][0][t] + aux_values_N[m + 1][k + 1][0][t]) / 2) * dx_2 -
								((aux_values_N[m][k + 1][0][t] + aux_values_N[m + 1][k + 1][0][t]) / 2) * dx_3 -
								((aux_values_N[m][k][0][t] + aux_values_N[m][k + 1][0][t]) / 2) * dx_4);
	}
	//------Second term
	if (smart_wds) {
		if ((((CdeP_vrt_1 + CdeP_vrt_2) > CdP_tsh * 0.5 * rho_inf * V_inf * V_inf * chord) || ((CdeP_vrt_1 + CdeP_vrt_2) < -CdP_tsh * 0.5 * rho_inf * V_inf * V_inf * chord)) && flags[m][k][2][t]) {
			AF_elem[m][k][17][t] = CdeP_vrt_2;
		} else if (flags[m][k][1][t]) {
			AF_elem[m][k][18][t] = CdeP_vrt_2;
		} else {
			AF_elem[m][k][34][t] = CdeP_vrt_2;
		}
	} else {
		if (flags[m][k][2][t]) {
			AF_elem[m][k][17][t] = CdeP_vrt_2;
		} else if (flags[m][k][1][t]) {
			AF_elem[m][k][18][t] = CdeP_vrt_2;
		} else {
			AF_elem[m][k][34][t] = CdeP_vrt_2;
		}
	}
	if (flux_scheme) {
		if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
			AF_elem[m][k][19][t] = 0;
		} else if (k != 0) {
			AF_elem[m][k][19][t]= (
				(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[m][k - 1][2][t] + (1 - a_1) * aux_values[m][k][2][t]) +
				(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][2][t] + (1 - a_2) * aux_values[m][k][2][t]) +
				(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][2][t] + (1 - a_3) * aux_values[m][k][2][t]) +
				(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][2][t] + (1 - a_4) * aux_values[m][k][2][t]));
		} else {
			if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
					AF_elem[m][k][19][t] = (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k][2][t] + (1 - a_1) * aux_values[m][k][2][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][2][t] + (1 - a_2) * aux_values[m][k][2][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][2][t] + (1 - a_3) * aux_values[m][k][2][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][2][t] + (1 - a_4) * aux_values[m][k][2][t]));
				} else {
					AF_elem[m][k][19][t] = (
						(rx_1 * dz_1 - rz_1 * dx_1) * (-a_1 * aux_values[m][k + 1][2][t] + (1 + a_1) * aux_values[m][k][2][t]) +
						(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][2][t] + (1 - a_2) * aux_values[m][k][2][t]) +
						(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][2][t] + (1 - a_3) * aux_values[m][k][2][t]) +
						(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][2][t] + (1 - a_4) * aux_values[m][k][2][t]));
				}
			} else {
				AF_elem[m][k][19][t] = (
					(rx_1 * dz_1 - rz_1 * dx_1) * (a_1 * aux_values[i_TOT - 2 - m][k][2][t] + (1 - a_1) * aux_values[m][k][2][t]) +
					(rx_2 * dz_2 - rz_2 * dx_2) * (a_2 * aux_values[m + 1][k][2][t] + (1 - a_2) * aux_values[m][k][2][t]) +
					(rx_3 * dz_3 - rz_3 * dx_3) * (a_3 * aux_values[m][k + 1][2][t] + (1 - a_3) * aux_values[m][k][2][t]) +
					(rx_4 * dz_4 - rz_4 * dx_4) * (a_4 * aux_values[m - 1][k][2][t] + (1 - a_4) * aux_values[m][k][2][t]));
			}
		}
	} else {
		AF_elem[m][k][19][t] = (
			(rx_1 * dz_1 - rz_1 * dx_1) * ((aux_values_N[m][k][2][t] + aux_values_N[m + 1][k][2][t]) / 2) +
			(rx_2 * dz_2 - rz_2 * dx_2) * ((aux_values_N[m + 1][k][2][t] + aux_values_N[m + 1][k + 1][2][t]) / 2) +
			(rx_3 * dz_3 - rz_3 * dx_3) * ((aux_values_N[m][k + 1][2][t] + aux_values_N[m + 1][k + 1][2][t]) / 2) +
			(rx_4 * dz_4 - rz_4 * dx_4) * ((aux_values_N[m][k][2][t] + aux_values_N[m][k + 1][2][t]) / 2));
	}
}

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
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
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
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
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
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
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
				if (( NS && (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))))) {
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

void compute_AF(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values,
				vector < vector < vector < vector < double > > > > &aux_values, vector < vector < vector < vector < double > > > > &aux_values_N,
				vector < vector < vector < vector < bool > > > > &flags, vector < vector < vector < vector < double > > > > &AF_elem,
				int t, int mD, int kD, int numVars, int compVarsN, int numAF, int numAF_ADD, int i_TOT, int64_t* kMax, bool NS, bool flux_scheme, int grad_scheme, bool ParForm, bool limit_integr,
				vector < bool > &near_ID, double rho_inf, double V_inf, double chord, int rho_id, int vx_id, int vz_id, vector < double > &X_POLE, vector < double > &Z_POLE,
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
							//To be implemented
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
							//To be implemented
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
							//To be implemented
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
					//To be implemented
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
					//To be implemented
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
		//SURFACE integrals over BODY boundary (only computed in case of eff. inviscid flow or moving-grid and inertial formulation)
		if ( (k==0) && ( ParForm || ((!NS)||(MOV_GRD && (U_FORM==0))) ) ) {
			for (int m = mD; m <i_TOT-1-mD; m++) {
				if ((near_ID[m]) || (near_ID[m + 1])) {
					cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
								rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
					//Vortex force additional term on body surface
					AERO[numAF+8][t]=AERO[numAF+8][t]-((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
															((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2)*dz_1;
					AERO[numAF+9][t]=AERO[numAF+9][t]+((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*
															((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2)*dx_1;
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
							//double r_prime_x = rx_1 - x_CoR+ampl_x * sin(freq_x*(time[t]-t0_x));
							//double r_prime_z = rz_1 - z_CoR+ampl_z * sin(freq_z*(time[t]-t0_z));
							double vx=NAN, vz=NAN;
							//To be corrected!!! rz and rx should be replaced by relative positions (x' vector components)
							vx = ((aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t])/2) - freq_x*freq_x*ampl_x*sin(freq_x*(time[t]-t0_x)) -
								freq_a*freq_a*ampl_a*sin(freq_a*(time[t]-t0_a))*rz_1 - pow(freq_a*ampl_a*cos(freq_a*(time[t]-t0_a)),2)*rx_1;
							vz = ((aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t])/2) - freq_z*freq_z*ampl_z*sin(freq_z*(time[t]-t0_z)) +
								freq_a*freq_a*ampl_a*sin(freq_a*(time[t]-t0_a))*rx_1 - pow(freq_a*ampl_a*cos(freq_a*(time[t]-t0_a)),2)*rz_1;
							AERO[numAF+32][t]=AERO[numAF+32][t]+rz_1*(dz_1*vz + dx_1*vx);
							AERO[numAF+33][t]=AERO[numAF+33][t]-rx_1*(dz_1*vz + dx_1*vx);
						}
					}
					//Unsteady body surface term (Formulation C)
					if (MOV_GRD && (U_FORM==0)) {	//Not computed in cases it should be exactely zero (to avoid inaccuracies due to non-zero relative velocity at trailing edge points)
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

int main(int argc, char** argv) {
	const double R_air = 287.058;
	const double gamma_air = 1.4;
	const double EPS = 1e-10;		//Tolerance for nodes matching in the 'wake line' (used to identify grid nodes on the body surface)
	const int IN_MARGIN=2;			//Number of cells added around selected integration domain to limit post-processed data in case of both m_SEL and k_SEL greater than zero.
	//Array cardinalities
	const int compVarsN = 2;
	const int compVarsC = 21;
	const int numFlags = 4;
	const int numAuxC = 16;
	const int numAuxN= 36;
	const int numAF = 41;
	const int numAF_ADD = 44;
	int numDBG=4;
	//Variables to be read from a case input file
	float m_SEL = 10;
	float k_SEL = 192;
	int TEC_OUT = 0;
	int vx_id = 6;
	int vz_id = 7;
	int rho_id = 4;
	int p_id = 4;
	int T_id = 5;
	int mut_id = 8;
	int k_id = 2;	//In case of turb.kin. energy not provided (as for 1-eq. turb. models), this default setting allows using y coordinates (zero) as turb. kin. energy.
	int flag_id = -1;
	int vgx_id = 2;	//In case of grid velocities not provided (as also in case of NO moving grid), this default setting allows using y coordinates (zero) as both x and z components of grid velocity.
	int vgz_id = 2;	//In case of grid velocities not provided (as also in case of NO moving grid), this default setting allows using y coordinates (zero) as both x and z components of grid velocity.
	double AoA = 0;
	bool NS=1;
	bool INC=0;
	bool RHO_VAR=0;
	double rho_inf = 1.225;
	double T_inf = 288.15;
	double V_inf = 100;
	double chord = 1.0;
	int pole_mode = 1;
	double X_POLE_IN = 0;
	double Z_POLE_IN = 0;
	double Fshock_tsh = 1.0;
	double Fshock_tail = 1.052;	//Shock function limit for viscous region cells to be considered in wave drag domain (has effect at shock tail)
	double K_bl = 1.1;
	double omega_tsh = 10;
	double entr_tsh = 0.001;
	//const double k_tsh=0.15;
	int shock_margin = 3;
	int shock_margin_drag_kmax = 2;
	int shock_margin_drag_mmin = 2;
	int shock_margin_drag_mmax = 2;
	//const double shock_top_margin=0.15;
	double shock_xds_margin = 0.0;
	//const double shock_xus_margin=0.01;
	int wall_margin = 20;
	int bl_margin_NF = 2;
	int bl_margin_FF_UP = 2;
	int bl_margin_FF_LOW = 2;
	double BL_sens_inf = 0;
	int grad_scheme = 0;
	bool flux_scheme = 1;
	bool limit_integr = 0;
	bool dom_decomp = 0;
	bool thd_methods = 0;
	bool ONERA_corr = 1;
	bool ParForm = 0;
	int mrho_form = 2;
	int Lamb_form = 0;
	bool smart_wds = 0;
	double CdP_tsh = 2.5e-6;	//Tolerance for elementary parassite wave drag coefficient of individual grid cells
	//Integration domain upper and lower bounds
	int k_LB=0;
	int m_LB=1;
	int k_UB=0;
	int m_UB=0;
	bool rescaled=0;
	//Unsteady variables
	bool MOV_GRD=0;		//This flag is necessary to distinguish whether the absolute velocity includes grid velocities or just represents the fluid velocity.
	int U_FORM=0;		// when (U_FORM==0) the absolute velocity is used; when (U_FORM==1 || U_FORM==2) the relative velocity is used;
	bool FT_FORM=0;
	double x_CoR=0;		double z_CoR=0;
	double ampl_x=0;	double freq_x=0;	double t0_x=0;
	double ampl_z=0;	double freq_z=0;	double t0_z=0;
	double ampl_a=0;	double freq_a=0;	double t0_a=0;
	int Nt=1;
	int every=1;
	//Chrono variables
	auto start = chrono::steady_clock::now();
	auto end = chrono::steady_clock::now();
	//other variables
	double p_inf = 0;
	double M_inf = 0;
	bool debug_mode = 0;

	//Set debug mode FLAG
	if (argc <= 3) {
		debug_mode = false;
	}
	else if (atof(argv[3]) == 1) {
			debug_mode = true;
	}

	//Read case and settings file
	ifstream case_file(argv[2]);
	if (case_file.is_open()) {
		string line;
		int i = 0;
		int m = 0;
		while (getline(case_file, line)) {
			i = line.find_last_of("\t ");
			m = line.find_first_of("\t ");
			if (line.substr(0, m) == "AoA") {
				AoA = atof(line.substr(i + 1).data());
			}
			if (line.substr(0, m) == "NS") {
				NS = atof(line.substr(i + 1).data());
			}
			if (line.substr(0, m) == "INC") {
				INC = atof(line.substr(i + 1).data());
			}
			if (line.substr(0, m) == "RHO_VAR") {
				RHO_VAR = atof(line.substr(i + 1).data());
			}
			if (line.substr(0, m) == "MOVGRD") {
				MOV_GRD = atof(line.substr(i + 1).data());
			}
			if (line.substr(0, m) == "UFORM") {
				U_FORM = atof(line.substr(i + 1).data());
			}
			if (line.substr(0, m) == "FT_FORM") {
				FT_FORM = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "VX_ID") {
				vx_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "VZ_ID") {
				vz_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "RHO_ID") {
				rho_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "T_ID") {
				T_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "P_ID") {
				p_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "MUT_ID") {
				mut_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "K_ID") {
				k_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "FLAG_ID") {
				flag_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "VGX_ID") {
				vgx_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "VGZ_ID") {
				vgz_id = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "V_INF") {
				V_inf = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "RHO_INF") {
				rho_inf = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "T_INF") {
				T_inf = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "CHORD") {
				chord = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "F_SHOCK_TSH") {
				Fshock_tsh = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "F_SHOCK_TAIL") {
				Fshock_tail = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "K_BL") {
				K_bl = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "OMEGA_TSH") {
				omega_tsh = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "ENTR_TSH") {
				entr_tsh = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "SHOCK_MARGIN") {
				shock_margin = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "SHOCK_MARGIN_DRAG_TOP") {
				shock_margin_drag_kmax = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "SHOCK_MARGIN_DRAG_MIN") {
				shock_margin_drag_mmin = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "SHOCK_MARGIN_DRAG_MAX") {
				shock_margin_drag_mmax = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "SHOCK_XDS_MARGIN") {
				shock_xds_margin = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "WALL_MARGIN") {
				wall_margin = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "BL_MARGIN_NF") {
				bl_margin_NF = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "BL_MARGIN_FF_UP") {
				bl_margin_FF_UP = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "BL_MARGIN_FF_LOW") {
				bl_margin_FF_LOW = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "CD_P_TSH") {
				CdP_tsh = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "I_SEL") {
				m_SEL = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "K_SEL") {
				k_SEL = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "GRAD_SCHEME") {
				grad_scheme = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "FLUX_SCHEME") {
				flux_scheme = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "ONERA_CORR") {
				ONERA_corr = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "PAR_FORM") {
				ParForm = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "LAMB_FORM") {
				Lamb_form= atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "LIMIT_INTEGR") {
				limit_integr = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "DOM_DECOMP") {
				dom_decomp = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "THD_METHODS") {
				thd_methods = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "MRHO_FORM") {
				mrho_form = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "SMART_WDS") {
				smart_wds = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "POLE_MODE") {
				pole_mode = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "X_POLE") {
				X_POLE_IN = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "Z_POLE") {
				Z_POLE_IN = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "TEC_OUT") {
				TEC_OUT = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "X_CoR") {
				x_CoR = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "Z_CoR") {
				z_CoR = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "AMPL_X") {
				ampl_x = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "FREQ_X") {
				freq_x = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "T0_X") {
				t0_x = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "AMPL_Z") {
				ampl_z = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "FREQ_Z") {
				freq_z = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "T0_Z") {
				t0_z = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "AMPL_A") {
				ampl_a = atof(line.substr(i + 1).data());
				ampl_a = ampl_a * M_PI / 180;
			}
			else if (line.substr(0, m) == "FREQ_A") {
				freq_a = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "T0_A") {
				t0_a = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "NT") {
				Nt = atof(line.substr(i + 1).data());
			}
			else if (line.substr(0, m) == "EVERY") {
				every = atof(line.substr(i + 1).data());
			}
		}
		case_file.close();
	}
	Nt=floor((Nt-1)/every)+1;
	p_inf= rho_inf * R_air * T_inf;
	M_inf = V_inf / (sqrt(gamma_air * R_air * T_inf));

	//Read Tecplot SZL files
	//assumes a C-type Multi-block structured mesh with possibilty of multiple blocks along the C-direction
	//C-direction index is assumed to coincide with the I index
	//All tecplot solution files, at different time steps, must have equal structure.
	start = chrono::steady_clock::now();
	void* inputFileHandle = NULL;
	int R;
	string line=argv[1];
	string file=line.data();
	int us = line.find_last_of("_");
	int pt = line.find_last_of(".");
	int32_t numZones = 0;
	if (debug_mode) {
		cout << "Reading Tecplot file: " << file.data() << "\n";
	}
	R = tecFileReaderOpen((char const*)(file.data()), &inputFileHandle);
	R = tecDataSetGetNumZones(inputFileHandle, &numZones);
	if (debug_mode) {
		cout << "-----ZONEs number is: " << numZones << "\n";
	}
	int64_t iMax[numZones], kMax[numZones];
		for (int i = 0; i < numZones; i++) {	//Initialization loop
			iMax[i] = 0;
			kMax[i] = 0;
		}
	int64_t tmp_max = 0;
	for (int i = 0; i < numZones; i++) {
		R = tecZoneGetIJK(inputFileHandle, i + 1, &iMax[i], &tmp_max, &kMax[i]);
		if (kMax[i]==1) {
			R = tecZoneGetIJK(inputFileHandle, i + 1, &iMax[i], &kMax[i], &tmp_max);
		}
	}
	if (debug_mode) {
		for (int i = 0; i < numZones; i++) {
			cout << "-------------iMax is: " << iMax[i] << "\n";
			cout << "-------------kMax is: " << kMax[i] << "\n";
		}
	}
	int32_t numVars = 0;
	R = tecDataSetGetNumVars(inputFileHandle, &numVars);
	int64_t numValues[numZones][numVars];
		for (int i = 0; i < numZones; i++) {	//Initialization loop
			for (int j = 0; j < numVars; j++) {
				numValues[i][j] = 0;
			}
		}
	string in_var_names[numVars];
		for (int i = 0; i < numVars; i++) {	//Initialization loop
			in_var_names[i] = "";
		}
	for (int j = 0; j < numVars; j++) {
		char* name = NULL;
		R = tecVarGetName(inputFileHandle, j + 1, &name);
		in_var_names[j] = name;
	}
	int64_t max_nv = 0;
	for (int i = 0; i < numZones; i++) {
		for (int j = 0; j < numVars; j++) {
			R = tecZoneVarGetNumValues(inputFileHandle, i + 1, j + 1, &numValues[i][j]);
			if (numValues[i][j] > max_nv) {
				max_nv = numValues[i][j];
			}
		}
	}
	if (debug_mode) {
		for (int i = 0; i < numZones; i++) {
			for (int j = 0; j < numVars; j++) {
				cout << "-----Values of zone " << i << " for variable " << j << " is:   " << numValues[i][j] << "\n";
			}
		}
		cout << "-----Maximum number of variable values is:   " << max_nv << "\n";
	}
	vector < vector < vector < double > > > values(numZones, vector < vector < double > > (numVars, vector < double > (max_nv,NAN)));
	vector < double > time(Nt,NAN);
	for (int i = 0; i < numZones; i++) {
		for (int j = 0; j < numVars; j++) {
			R = tecZoneVarGetDoubleValues(inputFileHandle, i + 1, j + 1, 1, numValues[i][j], &values[i][j][0]);
		}
	}
	R = tecZoneGetSolutionTime(inputFileHandle, 1, &time[0]);
	R = tecFileReaderClose(&inputFileHandle);
	//Reshape data based on a i-k-var 3D matrix. All blocks are merged into one
	int64_t i_TOT = iMax[0];
	for (int i = 1; i < numZones; i++) {
		i_TOT = i_TOT + iMax[i] - 1;
	}
	if (debug_mode) {
		cout << "-----Maximum number of global I index is:   " << i_TOT << "\n";
	}
	//Conditionally limit m and k cardinalities to save memory and computational time
	if ((numZones==1)&&((m_SEL>0)&&(k_SEL>0))) {
		i_TOT = (i_TOT - 2 * m_SEL) + 2 * IN_MARGIN;
		kMax[0] = k_SEL + 1 + IN_MARGIN;
		rescaled = true;
		if (debug_mode) {
			cout << "-----RESCALED DIMENSIONS TO SAVE RAM and CPU (INTEGRATION DOMAIN MARGIN = " << IN_MARGIN << ")\n";
			cout << "-----Maximum number of global I index is:   " << i_TOT << "\n";
			cout << "-----Maximum number k index is:             " << kMax[0] << "\n";
		}
	}
	//---Dynamically allocate and fill-in nodal values data structure
	vector < vector < vector < vector < double > > > > \
		nodal_values(i_TOT, vector < vector < vector < double > > > (kMax[0], vector < vector < double > > (numVars + compVarsN, vector < double > (Nt,0))));
	for (int j = 0; j < numVars; j++) {
		int riemp_i = 1;
		for (int i = 0; i < numZones; i++) {
			for (int k = 0; k < kMax[i]; k++) {
				if (rescaled) {
					for (int m = m_SEL-IN_MARGIN; m < iMax[0]-m_SEL+IN_MARGIN; m++) {
						nodal_values[m-(m_SEL-IN_MARGIN)][k][j][0] = values[i][j][k * iMax[0] + m];
					}
				} else {
					for (int m = 0; m < iMax[i]; m++) {
						nodal_values[riemp_i - 1 + m][k][j][0] = values[i][j][k * iMax[i] + m];
					}
				}
			}
			if (rescaled) {
				riemp_i = iMax[0]-2*m_SEL+2*IN_MARGIN;
			} else {
				riemp_i = riemp_i - 1 + iMax[i];
			}
			if (debug_mode) {
				cout << "-----Actual riemp-i is:   " << riemp_i << "\n";
			}
		}
	}

	//Successive time steps
	for (int t=1; t<Nt; t++) {
		stringstream middle;		
		middle << std::setw(pt-us-1) << std::setfill('0') << atof(line.substr(us+1, pt-us-1).data())+t*every;
		file=line.substr(0, us+1)+middle.str()+line.substr(pt);
		if (debug_mode) {
			cout << "Reading Tecplot file: " << file.data() << "\n";
		}
		R = tecFileReaderOpen((char const*)(file.data()), &inputFileHandle);		
		for (int i = 0; i < numZones; i++) {
			for (int j = 0; j < numVars; j++) {
				R = tecZoneVarGetDoubleValues(inputFileHandle, i + 1, j + 1, 1, numValues[i][j], &values[i][j][0]);
			}
		}
		R = tecZoneGetSolutionTime(inputFileHandle, 1, &time[t]);
		R = tecFileReaderClose(&inputFileHandle);
		for (int j = 0; j < numVars; j++) {
			int riemp_i = 1;
			for (int i = 0; i < numZones; i++) {
				for (int k = 0; k < kMax[i]; k++) {
					if (rescaled) {
						for (int m = m_SEL-IN_MARGIN; m < iMax[0]-m_SEL+IN_MARGIN; m++) {
							nodal_values[m-(m_SEL-IN_MARGIN)][k][j][t] = values[i][j][k * iMax[0] + m];
						}
					} else {
						for (int m = 0; m < iMax[i]; m++) {
							nodal_values[riemp_i - 1 + m][k][j][t] = values[i][j][k * iMax[i] + m];
						}
					}
				}
				if (rescaled) {
					riemp_i = iMax[0]-2*m_SEL+2*IN_MARGIN;
				} else {
					riemp_i = riemp_i - 1 + iMax[i];
				}
				if (debug_mode) {
					cout << "-----Actual riemp-i is:   " << riemp_i << "\n";
				}
			}
		}
	}
	values.clear();
	//Conditionally rescale m_SEL in case post-processed data were limited to save RAM and CPU.
	if (rescaled) {
		m_SEL=IN_MARGIN;
		if (debug_mode) {
			cout << "----------------------------------- " << "\n";
			cout << "-----Input m_SEL RESCALED to " << m_SEL << "\n";
			cout << "----------------------------------- " << "\n";
		}
	}
	//values.swap(vector < vector < vector < double > > >());
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of reading Tecplot files.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Conditionally replace absolute with relative velocities as working variables
	if ((U_FORM==1)||(U_FORM==2)) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0]; k++) {
				for (int m = 0; m < i_TOT; m++) {
					nodal_values[m][k][vx_id-1][t] = nodal_values[m][k][vx_id-1][t] - nodal_values[m][k][vgx_id-1][t];
					nodal_values[m][k][vz_id-1][t] = nodal_values[m][k][vz_id-1][t] - nodal_values[m][k][vgz_id-1][t];
				}
			}
		}
	}

	//---Data structures with dynamic allocation
	start = chrono::steady_clock::now();
	vector < vector < vector < vector < double > > > > \
		cell_values(i_TOT-1, vector < vector < vector < double > > > (kMax[0]-1, vector < vector < double > > (numVars + compVarsN + compVarsC, vector < double > (Nt,0))));
	vector < vector < vector < vector < bool > > > > \
		flags(i_TOT-1, vector < vector < vector < bool > > > (kMax[0]-1, vector < vector < bool > > (numFlags, vector < bool > (Nt,false))));
	vector < vector < vector < bool > > > flag_org(i_TOT-1, vector < vector < bool > > (kMax[0]-1, vector < bool > (Nt,false)));
	vector < vector < vector < vector < double > > > > \
		AF_elem(i_TOT-1, vector < vector < vector < double > > > (kMax[0]-1, vector < vector < double > > (numAF, vector < double > (Nt,0))));
	vector < vector < vector < vector < double > > > > \
		aux_values(i_TOT-1, vector < vector < vector < double > > > (kMax[0]-1, vector < vector < double > > (numAuxC, vector < double > (Nt,0))));
	vector < vector < vector < vector < double > > > > \
		aux_values_N(i_TOT, vector < vector < vector < double > > > (kMax[0], vector < vector < double > > (numAuxN, vector < double > (Nt,0))));

	double dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, a_1, a_2, a_3, a_4;
		dx_1 = 0; dz_1 = 0; dx_2 = 0; dz_2 = 0; dx_3 = 0; dz_3 = 0; dx_4 = 0; dz_4 = 0;
		rx_1 = 0; rz_1 = 0; rx_2 = 0; rz_2 = 0; rx_3 = 0; rz_3 = 0; rx_4 = 0; rz_4 = 0; a_1 = 0; a_2 = 0; a_3 = 0; a_4 = 0;
	vector < vector < double > > WLS(2,vector < double > (4,0));
	vector <bool> near_ID(i_TOT,false);
	vector < double > X_POLE(Nt,NAN);
	vector < double > Z_POLE(Nt,NAN);
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of allocating data structures.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Identify nodes on the body surface
	R=0;
	for (int m=0; m<i_TOT; m++) {
		if (((abs(nodal_values[m][0][0][0]-nodal_values[i_TOT-m-1][0][0][0])>EPS)||(abs(nodal_values[m][0][2][0]-nodal_values[i_TOT-m-1][0][2][0])>EPS))||(m==i_TOT-m-1)) {
			near_ID[m]=true;
			R++;
		}
	}
	if (debug_mode) {
		cout << "Identified " << R+1 << " cells on the body surface.\n";
	}

	//Setup integration domain upper and lower bounds
	start = chrono::steady_clock::now();
	int m_te=0;
	while ((m_te<=int(i_TOT/2))&&(!near_ID[m_te])) {
		m_te++;
	}
	m_te--;
	if (debug_mode) {
		cout << "-----m_te is: " << m_te << "\n";
	}
	//Find grid lines indices in case of domain selected in chord units away from the body
		if ((m_SEL<0)&&(k_SEL<0)) {
			double dist_old=1000000, dist;
			float m_SEL_tmp=-m_SEL;
			float k_SEL_tmp=-k_SEL;			
			// for (int k=0; k<kMax[0]; k++) {
			// 	for (int i=0; i<i_TOT; i++) {
			// 		dist=sqrt(pow(nodal_values[i][k][0][0]-(m_SEL_tmp+nodal_values[m_te][0][0][0]),2)+pow(nodal_values[i][k][2][0]+(k_SEL_tmp+nodal_values[m_te][0][2][0]),2));
			// 		if (dist<dist_old) {
			// 			m_SEL=i;
			// 			k_SEL=k;
			// 			dist_old=dist;
			// 		}
			// 	}
			// }
			dist_old=1000000;
			for (int i=0; i<=m_te; i++) {
				dist=abs(nodal_values[i][0][0][0]-(m_SEL_tmp+nodal_values[m_te][0][0][0]));
				if (dist<dist_old) {
					m_SEL=i;
					dist_old=dist;
				}
			}
			dist_old=1000000;
			for (int k=0; k<kMax[0]-1; k++) {
				dist=abs(nodal_values[m_SEL][k][2][0]+(k_SEL_tmp+nodal_values[m_te][0][2][0]));
				if (dist<dist_old) {
					k_SEL=k;
					dist_old=dist;
				}
			}
			if (m_SEL<1) {
				m_SEL=1;
			}
			if (k_SEL<1) {
				k_SEL=1;
			}
			if (debug_mode) {
				cout << "-----Selected domain based on provided distance values:   m_SEL=" << m_SEL << " ; k_SEL=" << k_SEL << "\n";
			}
		}
	if (m_SEL>0) {
		m_LB=int(m_SEL);
		m_UB=int(m_SEL)+1;
	} else {
		m_UB=i_TOT-1;
	}
	if (k_SEL>0) {
		k_LB=int(k_SEL);
		k_UB=int(k_SEL)+1;
	} else {
		k_UB=kMax[0]-1;
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of setup integration domain upper and lower bounds.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//----------------------------------------POST-PROCESSING----------------------------------------------------
	if (INC) {
		rho_id=p_id;
	}
	//---Solutions rotation (around grid origin) to align with wind axes
	start = chrono::steady_clock::now();
	if (AoA!=0) {	//To save computational time on zero AoA cases.
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0]; k++) {
				for (int m = 0; m < i_TOT; m++) {
					double x_tmp = nodal_values[m][k][0][t];
					double z_tmp = nodal_values[m][k][2][t];
					nodal_values[m][k][0][t] = x_tmp * cos(AoA * M_PI / 180) + z_tmp * sin(AoA * M_PI / 180);
					nodal_values[m][k][2][t] = -x_tmp * sin(AoA * M_PI / 180) + z_tmp * cos(AoA * M_PI / 180);
					x_tmp = nodal_values[m][k][vx_id - 1][t];
					z_tmp = nodal_values[m][k][vz_id - 1][t];
					nodal_values[m][k][vx_id - 1][t] = x_tmp * cos(AoA * M_PI / 180) + z_tmp * sin(AoA * M_PI / 180);
					nodal_values[m][k][vz_id - 1][t] = -x_tmp * sin(AoA * M_PI / 180) + z_tmp * cos(AoA * M_PI / 180);
				}
			}
		}
	}
	//---Set pole location and conditionally rotate it.
	switch (pole_mode) {
		case 1:		//Pole location provided in original inertial grid coordinates: need to be rotated
			for (int t=0; t<Nt; t++) {
				X_POLE[t] = X_POLE_IN * cos(AoA * M_PI / 180) + Z_POLE_IN * sin(AoA * M_PI / 180);
				Z_POLE[t] = -X_POLE_IN * sin(AoA * M_PI / 180) + Z_POLE_IN * cos(AoA * M_PI / 180);
			}
			cout << "---Pole fixed at (" << X_POLE[0] << "," << Z_POLE[0] << ") in post-processor inertial coordinates.\n";
			break;
		case 2:		//Pole location provided in post-processor inertial coordinates: take as they are
			for (int t=0; t<Nt; t++) {
				X_POLE[t] = X_POLE_IN;
				Z_POLE[t] = Z_POLE_IN;
			}
			cout << "---Pole fixed at (" << X_POLE[0] << "," << Z_POLE[0] << ") in post-processor inertial coordinates.\n";
			break;
		case 3:		//Pole location provided as m and k grid lines indices (starting from 0). This is a non-inertial pole definition mode.
			for (int t=0; t<Nt; t++) {
				X_POLE[t] = nodal_values[int(X_POLE_IN)][int(Z_POLE_IN)][0][t];
				Z_POLE[t] = nodal_values[int(X_POLE_IN)][int(Z_POLE_IN)][2][t];
			}
			cout << "---Pole initially set at (" << X_POLE[0] << "," << Z_POLE[0] << ") in post-processor inertial coordinates.\n";
			break;
		case 4:		//Pole location provided in post-processor non-inertial coordinates: apply input grid motion to get inertial pole coordinates
			for (int t=0; t<Nt; t++) {
				X_POLE[t] = x_CoR + ampl_x*sin(freq_x*(time[t]-t0_x)) + X_POLE_IN*cos(ampl_a*sin(freq_a*(time[t]-t0_a))) + Z_POLE_IN*sin(ampl_a*sin(freq_a*(time[t]-t0_a)));
				Z_POLE[t] = z_CoR + ampl_z*sin(freq_z*(time[t]-t0_z)) - X_POLE_IN*sin(ampl_a*sin(freq_a*(time[t]-t0_a))) + Z_POLE_IN*cos(ampl_a*sin(freq_a*(time[t]-t0_a)));
			}
			cout << "-----Pole initially set at (" << X_POLE[0] << "," << Z_POLE[0] << ") in post-processor inertial coordinates.\n";
			break;
		case 0:		//Pole location is automatically set on the wake plane of the integraion domain. Input data on X_POLE_IN and Z_POLE_IN are ignored.
			//This case is not implemented for the unsteady post-processor.
			break;
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of solutions and pole rotation.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Compute additional variables at GRID NODES
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		for (int k = 0; k < kMax[0]; k++) {
			for (int m = 0; m < i_TOT; m++) {
				//Static pressure
					if (INC) {
						nodal_values[m][k][numVars][t] = nodal_values[m][k][p_id - 1][t];
						if (RHO_VAR) {
							nodal_values[m][k][rho_id - 1][t] = nodal_values[m][k][numVars][t] / ( R_air * nodal_values[m][k][T_id - 1][t] );
						} else {
							nodal_values[m][k][rho_id - 1][t] = rho_inf;
						}
					} else {
						nodal_values[m][k][numVars][t] = nodal_values[m][k][rho_id - 1][t] * R_air * nodal_values[m][k][T_id - 1][t];
					}
				//Kinetic energy per unit mass
				nodal_values[m][k][numVars + 1][t] = (pow(nodal_values[m][k][vx_id - 1][t], 2) + pow(nodal_values[m][k][vz_id - 1][t], 2)) / 2;
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of computing additional variables at grid nodes.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Average values of nodal input variables at cells centres
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		for (int j = 0; j < numVars + compVarsN; j++) {
			for (int k = 0; k < kMax[0] - 1; k++) {
				for (int m = 0; m < i_TOT - 1; m++) {
					if (j!=1) {	//Skip zero y-coordinates
						if (j == flag_id-1) {
							if (((nodal_values[m][k][j][t] == 0) && ((nodal_values[m+1][k][j][t] == 0) || (nodal_values[m+1][k+1][j][t] == 0) || (nodal_values[m][k+1][j][t] == 0))) ||
								((nodal_values[m+1][k][j][t] == 0) && ((nodal_values[m][k][j][t] == 0) || (nodal_values[m+1][k+1][j][t] == 0) || (nodal_values[m][k+1][j][t] == 0))) ||
								((nodal_values[m+1][k+1][j][t] == 0) && ((nodal_values[m][k][j][t] == 0) || (nodal_values[m+1][k][j][t] == 0) || (nodal_values[m][k+1][j][t] == 0))) ||
								((nodal_values[m][k+1][j][t] == 0) && ((nodal_values[m][k][j][t] == 0) || (nodal_values[m+1][k][j][t] == 0) || (nodal_values[m+1][k+1][j][t] == 0)))) {
									cell_values[m][k][j][t] = 0;
								} else if (((nodal_values[m][k][j][t] == 1) && ((nodal_values[m+1][k][j][t] == 1) || (nodal_values[m+1][k+1][j][t] == 1) || (nodal_values[m][k+1][j][t] == 1))) ||
										((nodal_values[m+1][k][j][t] == 1) && ((nodal_values[m][k][j][t] == 1) || (nodal_values[m+1][k+1][j][t] == 1) || (nodal_values[m][k+1][j][t] == 1))) ||
										((nodal_values[m+1][k+1][j][t] == 1) && ((nodal_values[m][k][j][t] == 1) || (nodal_values[m+1][k][j][t] == 1) || (nodal_values[m][k+1][j][t] == 1))) ||
										((nodal_values[m][k+1][j][t] == 1) && ((nodal_values[m][k][j][t] == 1) || (nodal_values[m+1][k][j][t] == 1) || (nodal_values[m+1][k+1][j][t] == 1)))) {
									cell_values[m][k][j][t] = 1;
								} else if (((nodal_values[m][k][j][t] + nodal_values[m + 1][k][j][t] + nodal_values[m][k + 1][j][t] + nodal_values[m + 1][k + 1][j][t]) / 4.0) >= 1.5) {
									cell_values[m][k][j][t] = 2;
								}
						} else {
							cell_values[m][k][j][t] = (nodal_values[m][k][j][t] + nodal_values[m + 1][k][j][t] + nodal_values[m][k + 1][j][t] + nodal_values[m + 1][k + 1][j][t]) / 4.0;
						}
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of averaging values of nodal input variables at cells centres.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Compute downstream grid line index of shock wake region
	int m_ds=0;
	while ((nodal_values[m_ds][0][0][0]-nodal_values[m_te][0][0][0])>shock_xds_margin) {
		m_ds++;
	}
	if (debug_mode) {
		cout << "-----m_ds is: " << m_ds << "\n";
	}

	//---Compute additional variables at CELLS CENTRES
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		for (int k = 0; k <kMax[0]-1; k++) {		//Cells area and centroid
			for (int m = 0; m <i_TOT-1; m++) {
				//Cell centroid coordinates
				//double x_FC = 0;
				//double z_FC = 0;
				double x_FC=(nodal_values[m][k][0][t]+nodal_values[m+1][k][0][t]+nodal_values[m+1][k+1][0][t]+nodal_values[m][k+1][0][t])/4;
				double z_FC=(nodal_values[m][k][2][t]+nodal_values[m+1][k][2][t]+nodal_values[m+1][k+1][2][t]+nodal_values[m][k+1][2][t])/4;
				// double x_T1 = (nodal_values[m][k][0] + nodal_values[m + 1][k][0] + nodal_values[m + 1][k + 1][0]) / 3;
				// double z_T1 = (nodal_values[m][k][2] + nodal_values[m + 1][k][2] + nodal_values[m + 1][k + 1][2]) / 3;
				// double x_T2 = (nodal_values[m][k][0] + nodal_values[m + 1][k + 1][0] + nodal_values[m][k + 1][0]) / 3;
				// double z_T2 = (nodal_values[m][k][2] + nodal_values[m + 1][k + 1][2] + nodal_values[m][k + 1][2]) / 3;
				// double x_T3 = (nodal_values[m][k][0] + nodal_values[m + 1][k][0] + nodal_values[m][k + 1][0]) / 3;
				// double z_T3 = (nodal_values[m][k][2] + nodal_values[m + 1][k][2] + nodal_values[m][k + 1][2]) / 3;
				// double x_T4 = (nodal_values[m + 1][k][0] + nodal_values[m + 1][k + 1][0] + nodal_values[m][k + 1][0]) / 3;
				// double z_T4 = (nodal_values[m + 1][k][2] + nodal_values[m + 1][k + 1][2] + nodal_values[m][k + 1][2]) / 3;
				// bool tmp = 0;
				// tmp = LineLineIntersect(x_T2, z_T2, x_T1, z_T1, x_T3, z_T3, x_T4, z_T4, x_FC, z_FC);
				// if (!tmp) {
				// 	cout << "ERROR\n";
				// }

				// typedef dense2D<double>  matrix_type;
				// matrix_type A_T12(2, 2), A_T34(2, 2), A_T(2, 2);
				// dense_vector<double> x_T12(2), b_T12(2), x_T34(2), b_T34(2), x_T(2), b_T(2);
				// //Line joining centres of triangles 1 and 2
				// 	A_T12(0,0)=x_T1;	A_T12(1,0)=x_T2;	A_T12(0,1)=1;	A_T12(1,1)=1;
				// 	pc::ilu_0<matrix_type> L_T12(A_T12);
				// 	//pc::identity<matrix_type> R_T12(A_T12);
				// 	b_T12(0)=z_T1; b_T12(1)=z_T2; x_T12(0)= 0; x_T12(1)= z_T1; 
				// 	noisy_iteration<double> iter_T12(A_T12*x_T12-b_T12, 50, 1.e-15);
				// 	bicgstab_2(A_T12, x_T12, b_T12, L_T12, iter_T12);
				// 	//qmr(A_T12, x_T12, b_T12, L_T12, R_T12, iter_T12);
				// //Line joining centres of triangles 3 and 4
				// 	A_T34(0,0)=x_T3;	A_T34(1,0)=x_T4;	A_T34(0,1)=1;	A_T34(1,1)=1;
				// 	pc::ilu_0<matrix_type> L_T34(A_T34);
				// 	//pc::identity<matrix_type> R_T34(A_T34);
				// 	b_T34(0)=z_T3; b_T34(1)=z_T4; x_T34(0)= 0; x_T34(1)= z_T3; 
				// 	noisy_iteration<double> iter_T34(A_T34*x_T34-b_T34, 50, 1.e-15);			
				// 	bicgstab_2(A_T34, x_T34, b_T34, L_T34, iter_T34);
				// 	//qmr(A_T34, x_T34, b_T34, L_T34, R_T34, iter_T34);
				// //Line-Line intersection to find quadrilateral cell centroid
				// 	A_T(0,0)=-x_T12(0);	A_T(1,0)=-x_T34(0);	A_T(0,1)=1;	A_T(1,1)=1;
				// 	pc::ilu_0<matrix_type> L_T(A_T);
				// 	//pc::identity<matrix_type> R_T(A_T);
				// 	b_T(0)=x_T12(1); b_T(1)=x_T34(1); x_T(0)=x_FC; x_T(1)=z_FC;
				// 	noisy_iteration<double> iter_T(A_T*x_T-b_T, 50, 1.e-15);
				// 	bicgstab_2(A_T, x_T, b_T, L_T, iter_T);
				// 	//qmr(A_T, x_T, b_T, L_T, R_T, iter_T);
				// x_FC=x_T(0);
				// z_FC=x_T(1);
							
				//Cell area (only computed at the first time step as the grid is not deforming)
				if (t==0) {
					cell_values[m][k][numVars+compVarsN][t]=
						S_tria(nodal_values[m][k][0][t], nodal_values[m][k][2][t], nodal_values[m+1][k][0][t], nodal_values[m+1][k][2][t], x_FC, z_FC)+
						S_tria(nodal_values[m+1][k][0][t], nodal_values[m+1][k][2][t], nodal_values[m+1][k+1][0][t], nodal_values[m+1][k+1][2][t], x_FC, z_FC)+
						S_tria(nodal_values[m+1][k+1][0][t], nodal_values[m+1][k+1][2][t], nodal_values[m][k+1][0][t], nodal_values[m][k+1][2][t], x_FC, z_FC)+
						S_tria(nodal_values[m][k+1][0][t], nodal_values[m][k+1][2][t], nodal_values[m][k][0][t], nodal_values[m][k][2][t], x_FC, z_FC);
						if (cell_values[m][k][numVars + compVarsN][t] < 0) {
							cout << "NEGATIVE area found at time step " << t << "and cell [" << m << "," << k << "] !   -   Value is: " << cell_values[m][k][numVars + compVarsN][t] << "\n";
						}
				} else {
					cell_values[m][k][numVars+compVarsN][t]=cell_values[m][k][numVars+compVarsN][0];
				}

				// //Cell area
				// cell_values[m][k][numVars+compVarsN][t]=0.5*((nodal_values[m][k][0][t]*nodal_values[m+1][k][2][t]-nodal_values[m+1][k][0][t]*nodal_values[m][k][2][t])+
				// 										  (nodal_values[m+1][k][0][t]*nodal_values[m+1][k+1][2][t]-nodal_values[m+1][k+1][0][t]*nodal_values[m+1][k][2][t])+
				// 										  (nodal_values[m+1][k+1][0][t]*nodal_values[m][k+1][2][t]-nodal_values[m][k+1][0][t]*nodal_values[m+1][k+1][2][t])+
				// 										  (nodal_values[m][k+1][0][t]*nodal_values[m][k][2][t]-nodal_values[m][k][0][t]*nodal_values[m][k+1][2][t]));
				// //Cell centroid coordinates
				// x_FC=(1/(6*cell_values[m][k][numVars+compVarsN][t]))*
				// 	 ((nodal_values[m][k][0][t]+nodal_values[m+1][k][0][t])*(nodal_values[m][k][0][t]*nodal_values[m+1][k][2][t]-nodal_values[m+1][k][0][t]*nodal_values[m][k][2][t])+
				// 	 (nodal_values[m+1][k][0][t]+nodal_values[m+1][k+1][0][t])*(nodal_values[m+1][k][0][t]*nodal_values[m+1][k+1][2][t]-nodal_values[m+1][k+1][0][t]*nodal_values[m+1][k][2][t])+
				// 	 (nodal_values[m+1][k+1][0][t]+nodal_values[m][k+1][0][t])*(nodal_values[m+1][k+1][0][t]*nodal_values[m][k+1][2][t]-nodal_values[m][k+1][0][t]*nodal_values[m+1][k+1][2][t])+
				// 	 (nodal_values[m][k+1][0][t]+nodal_values[m][k][0][t])*(nodal_values[m][k+1][0][t]*nodal_values[m][k][2][t]-nodal_values[m][k][0][t]*nodal_values[m][k+1][2][t]));
				// z_FC=(1/(6*cell_values[m][k][numVars+compVarsN][t]))*
				// 	 ((nodal_values[m][k][2][t]+nodal_values[m+1][k][2][t])*(nodal_values[m][k][0][t]*nodal_values[m+1][k][2][t]-nodal_values[m+1][k][0][t]*nodal_values[m][k][2][t])+
				// 	 (nodal_values[m+1][k][2][t]+nodal_values[m+1][k+1][2][t])*(nodal_values[m+1][k][0][t]*nodal_values[m+1][k+1][2][t]-nodal_values[m+1][k+1][0][t]*nodal_values[m+1][k][2][t])+
				// 	 (nodal_values[m+1][k+1][2][t]+nodal_values[m][k+1][2][t])*(nodal_values[m+1][k+1][0][t]*nodal_values[m][k+1][2][t]-nodal_values[m][k+1][0][t]*nodal_values[m+1][k+1][2][t])+
				// 	 (nodal_values[m][k+1][2][t]+nodal_values[m][k][2][t])*(nodal_values[m][k+1][0][t]*nodal_values[m][k][2][t]-nodal_values[m][k][0][t]*nodal_values[m][k+1][2][t]));
					
				//Cell centroid coordinates
				// double x_FC_old=x_FC, z_FC_old=z_FC;
				// x_FC=(1/cell_values[m][k][numVars+compVarsN])*(
				// 	S_tria(nodal_values[m][k][0], nodal_values[m][k][2], nodal_values[m+1][k][0], nodal_values[m+1][k][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m][k][0]+nodal_values[m+1][k][0]+x_FC_old)/3+
				// 	S_tria(nodal_values[m+1][k][0], nodal_values[m+1][k][2], nodal_values[m+1][k+1][0], nodal_values[m+1][k+1][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m+1][k][0]+nodal_values[m+1][k+1][0]+x_FC_old)/3+
				// 	S_tria(nodal_values[m+1][k+1][0], nodal_values[m+1][k+1][2], nodal_values[m][k+1][0], nodal_values[m][k+1][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m+1][k+1][0]+nodal_values[m][k+1][0]+x_FC_old)/3+
				// 	S_tria(nodal_values[m][k+1][0], nodal_values[m][k+1][2], nodal_values[m][k][0], nodal_values[m][k][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m][k+1][0]+nodal_values[m][k][0]+x_FC_old)/3);
				// z_FC=(1/cell_values[m][k][numVars+compVarsN])*(
				// 	S_tria(nodal_values[m][k][0], nodal_values[m][k][2], nodal_values[m+1][k][0], nodal_values[m+1][k][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m][k][2]+nodal_values[m+1][k][2]+z_FC_old)/3+
				// 	S_tria(nodal_values[m+1][k][0], nodal_values[m+1][k][2], nodal_values[m+1][k+1][0], nodal_values[m+1][k+1][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m+1][k][2]+nodal_values[m+1][k+1][2]+z_FC_old)/3+
				// 	S_tria(nodal_values[m+1][k+1][0], nodal_values[m+1][k+1][2], nodal_values[m][k+1][0], nodal_values[m][k+1][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m+1][k+1][2]+nodal_values[m][k+1][2]+z_FC_old)/3+
				// 	S_tria(nodal_values[m][k+1][0], nodal_values[m][k+1][2], nodal_values[m][k][0], nodal_values[m][k][2], x_FC_old, z_FC_old)*
				// 	(nodal_values[m][k+1][2]+nodal_values[m][k][2]+z_FC_old)/3);

				//Cells centroids coordinates
				cell_values[m][k][numVars+compVarsN+15][t]=(x_FC-X_POLE[t]);
				cell_values[m][k][numVars+compVarsN+16][t]=(z_FC-Z_POLE[t]);
			}
		}
		//---Spatial loops closed and re-opened to allow WLS gradient calculation using updated centroids coordinates)
		for (int k = 0; k <kMax[0]-1; k++) {		//Spatial gradients and additional variables
			for (int m = 0; m <i_TOT-1; m++) {
				cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
						rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
				//Dynamic (molecular) viscosity
				cell_values[m][k][numVars+compVarsN+1][t]=(0.000001458*sqrt(cell_values[m][k][T_id-1][t]))/(1+110.4/cell_values[m][k][T_id-1][t]);
				//Spatial gradients
				switch (grad_scheme) {
					case 0:		//Green-Gauss (nodes-based)
						//---Velocity
						cell_values[m][k][numVars+compVarsN+2][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							(((nodal_values[m+1][k][vx_id-1][t]+nodal_values[m][k][vx_id-1][t])/2)*dz_1+
							((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_2+
							((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
							((nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_4);
						cell_values[m][k][numVars+compVarsN+3][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							((-(nodal_values[m+1][k][vx_id-1][t]+nodal_values[m][k][vx_id-1][t])/2)*dx_1+
								(-(nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dx_2+
								(-(nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dx_3+
								(-(nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dx_4);
						cell_values[m][k][numVars+compVarsN+4][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							(((nodal_values[m+1][k][vz_id-1][t]+nodal_values[m][k][vz_id-1][t])/2)*dz_1+
								((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*dz_2+
								((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*dz_3+
								((nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*dz_4);
						cell_values[m][k][numVars+compVarsN+5][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							((-(nodal_values[m+1][k][vz_id-1][t]+nodal_values[m][k][vz_id-1][t])/2)*dx_1+
								(-(nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*dx_2+
								(-(nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*dx_3+
								(-(nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*dx_4);
						//---Pressure
						cell_values[m][k][numVars+compVarsN+6][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							(((nodal_values[m+1][k][numVars][t]+nodal_values[m][k][numVars][t])/2)*dz_1+
								((nodal_values[m+1][k+1][numVars][t]+nodal_values[m+1][k][numVars][t])/2)*dz_2+
								((nodal_values[m][k+1][numVars][t]+nodal_values[m+1][k+1][numVars][t])/2)*dz_3+
								((nodal_values[m][k][numVars][t]+nodal_values[m][k+1][numVars][t])/2)*dz_4);
						cell_values[m][k][numVars+compVarsN+7][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							((-(nodal_values[m+1][k][numVars][t]+nodal_values[m][k][numVars][t])/2)*dx_1+
								(-(nodal_values[m+1][k+1][numVars][t]+nodal_values[m+1][k][numVars][t])/2)*dx_2+
								(-(nodal_values[m][k+1][numVars][t]+nodal_values[m+1][k+1][numVars][t])/2)*dx_3+
								(-(nodal_values[m][k][numVars][t]+nodal_values[m][k+1][numVars][t])/2)*dx_4);
						//---Density
						cell_values[m][k][numVars+compVarsN+8][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)*dz_1+
								((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*dz_2+
								((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*dz_3+
								((nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*dz_4);
						cell_values[m][k][numVars+compVarsN+9][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							((-(nodal_values[m+1][k][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)*dx_1+
								(-(nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)*dx_2+
								(-(nodal_values[m][k+1][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)*dx_3+
								(-(nodal_values[m][k][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)*dx_4);
						//---Kinetic energy per unit mass
						cell_values[m][k][numVars+compVarsN+10][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							(((nodal_values[m+1][k][numVars+1][t]+nodal_values[m][k][numVars+1][t])/2)*dz_1+
								((nodal_values[m+1][k+1][numVars+1][t]+nodal_values[m+1][k][numVars+1][t])/2)*dz_2+
								((nodal_values[m][k+1][numVars+1][t]+nodal_values[m+1][k+1][numVars+1][t])/2)*dz_3+
								((nodal_values[m][k][numVars+1][t]+nodal_values[m][k+1][numVars+1][t])/2)*dz_4);
						cell_values[m][k][numVars+compVarsN+11][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							((-(nodal_values[m+1][k][numVars+1][t]+nodal_values[m][k][numVars+1][t])/2)*dx_1+
								(-(nodal_values[m+1][k+1][numVars+1][t]+nodal_values[m+1][k][numVars+1][t])/2)*dx_2+
								(-(nodal_values[m][k+1][numVars+1][t]+nodal_values[m+1][k+1][numVars+1][t])/2)*dx_3+
								(-(nodal_values[m][k][numVars+1][t]+nodal_values[m][k+1][numVars+1][t])/2)*dx_4);
						//---Turbulent kinetic energy per unit mass (assign it to auxiliary variables)
						aux_values[m][k][12][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							(((nodal_values[m+1][k][k_id-1][t]+nodal_values[m][k][k_id-1][t])/2)*dz_1+
								((nodal_values[m+1][k+1][k_id-1][t]+nodal_values[m+1][k][k_id-1][t])/2)*dz_2+
								((nodal_values[m][k+1][k_id-1][t]+nodal_values[m+1][k+1][k_id-1][t])/2)*dz_3+
								((nodal_values[m][k][k_id-1][t]+nodal_values[m][k+1][k_id-1][t])/2)*dz_4);
						aux_values[m][k][13][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							((-(nodal_values[m+1][k][k_id-1][t]+nodal_values[m][k][k_id-1][t])/2)*dx_1+
								(-(nodal_values[m+1][k+1][k_id-1][t]+nodal_values[m+1][k][k_id-1][t])/2)*dx_2+
								(-(nodal_values[m][k+1][k_id-1][t]+nodal_values[m+1][k+1][k_id-1][t])/2)*dx_3+
								(-(nodal_values[m][k][k_id-1][t]+nodal_values[m][k+1][k_id-1][t])/2)*dx_4);
						break;
					case 1:		//FD method (structured grids)
						//Mass density, fluid velocity,  specific kinetic energy and turbulent kinetic energy gradients are computed at grid nodes in the following sections, then averaged to cells centres.
						//Static pressure gradient is always computed using the G-G nodes-based scheme, as it is only used in the shock sensor.
						cell_values[m][k][numVars+compVarsN+6][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							(((nodal_values[m+1][k][numVars][t]+nodal_values[m][k][numVars][t])/2)*dz_1+
								((nodal_values[m+1][k+1][numVars][t]+nodal_values[m+1][k][numVars][t])/2)*dz_2+
								((nodal_values[m][k+1][numVars][t]+nodal_values[m+1][k+1][numVars][t])/2)*dz_3+
								((nodal_values[m][k][numVars][t]+nodal_values[m][k+1][numVars][t])/2)*dz_4);
						cell_values[m][k][numVars+compVarsN+7][t]=(1/cell_values[m][k][numVars+compVarsN][t])*
							((-(nodal_values[m+1][k][numVars][t]+nodal_values[m][k][numVars][t])/2)*dx_1+
								(-(nodal_values[m+1][k+1][numVars][t]+nodal_values[m+1][k][numVars][t])/2)*dx_2+
								(-(nodal_values[m][k+1][numVars][t]+nodal_values[m+1][k+1][numVars][t])/2)*dx_3+
								(-(nodal_values[m][k][numVars][t]+nodal_values[m][k+1][numVars][t])/2)*dx_4);
						break;
					case 2:		//Weighted Least Squares method
						double D_1, D_2, D_3, D_4;
						//---X-Velocity
							D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
							if (k != 0) {
								D_1=cell_values[m][k - 1][vx_id-1][t]-cell_values[m][k][vx_id-1][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									D_1=-2*cell_values[m][k][vx_id-1][t];
								} else {
									D_1=cell_values[i_TOT - 2 - m][k][vx_id-1][t]-cell_values[m][k][vx_id-1][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=cell_values[m + 1][k][vx_id-1][t]-cell_values[m][k][vx_id-1][t];
							}
							if (k != kMax[0] - 2) {
								D_3=cell_values[m][k + 1][vx_id-1][t]-cell_values[m][k][vx_id-1][t];
							}
							if (m != 0) {
								D_4=cell_values[m - 1][k][vx_id-1][t]-cell_values[m][k][vx_id-1][t];
							}
							cell_values[m][k][numVars+compVarsN+2][t]=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
							cell_values[m][k][numVars+compVarsN+3][t]=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						//---Z-Velocity
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=cell_values[m][k - 1][vz_id-1][t]-cell_values[m][k][vz_id-1][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									D_1=-2*cell_values[m][k][vz_id-1][t];
								} else {
									D_1=cell_values[i_TOT - 2 - m][k][vz_id-1][t]-cell_values[m][k][vz_id-1][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=cell_values[m + 1][k][vz_id-1][t]-cell_values[m][k][vz_id-1][t];
							}
							if (k != kMax[0] - 2) {
								D_3=cell_values[m][k + 1][vz_id-1][t]-cell_values[m][k][vz_id-1][t];
							}
							if (m != 0) {
								D_4=cell_values[m - 1][k][vz_id-1][t]-cell_values[m][k][vz_id-1][t];
							}
							cell_values[m][k][numVars+compVarsN+4][t]=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
							cell_values[m][k][numVars+compVarsN+5][t]=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						//---Pressure
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=cell_values[m][k - 1][numVars][t]-cell_values[m][k][numVars][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									//no action, D_1 equal zero (initialization) is ok
								} else {
									D_1=cell_values[i_TOT - 2 - m][k][numVars][t]-cell_values[m][k][numVars][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=cell_values[m + 1][k][numVars][t]-cell_values[m][k][numVars][t];
							}
							if (k != kMax[0] - 2) {
								D_3=cell_values[m][k + 1][numVars][t]-cell_values[m][k][numVars][t];
							}
							if (m != 0) {
								D_4=cell_values[m - 1][k][numVars][t]-cell_values[m][k][numVars][t];
							}
						cell_values[m][k][numVars+compVarsN+6][t]=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
						cell_values[m][k][numVars+compVarsN+7][t]=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						//---Density
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=cell_values[m][k - 1][rho_id-1][t]-cell_values[m][k][rho_id-1][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									D_1=-(cell_values[m][k + 1][rho_id-1][t]-cell_values[m][k][rho_id-1][t]);
								} else {
									D_1=cell_values[i_TOT - 2 - m][k][rho_id-1][t]-cell_values[m][k][rho_id-1][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=cell_values[m + 1][k][rho_id-1][t]-cell_values[m][k][rho_id-1][t];
							}
							if (k != kMax[0] - 2) {
								D_3=cell_values[m][k + 1][rho_id-1][t]-cell_values[m][k][rho_id-1][t];
							}
							if (m != 0) {
								D_4=cell_values[m - 1][k][rho_id-1][t]-cell_values[m][k][rho_id-1][t];
							}
						cell_values[m][k][numVars+compVarsN+8][t]=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
						cell_values[m][k][numVars+compVarsN+9][t]=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						//---Kinetic energy per unit mass
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=cell_values[m][k - 1][numVars+1][t]-cell_values[m][k][numVars+1][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									D_1=-2*cell_values[m][k][numVars+1][t];
								} else {
									D_1=cell_values[i_TOT - 2 - m][k][numVars+1][t]-cell_values[m][k][numVars+1][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=cell_values[m + 1][k][numVars+1][t]-cell_values[m][k][numVars+1][t];
							}
							if (k != kMax[0] - 2) {
								D_3=cell_values[m][k + 1][numVars+1][t]-cell_values[m][k][numVars+1][t];
							}
							if (m != 0) {
								D_4=cell_values[m - 1][k][numVars+1][t]-cell_values[m][k][numVars+1][t];
							}
						cell_values[m][k][numVars+compVarsN+10][t]=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
						cell_values[m][k][numVars+compVarsN+11][t]=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						//---Turbulent kinetic energy per unit mass
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=cell_values[m][k - 1][k_id-1][t]-cell_values[m][k][k_id-1][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									D_1=-(cell_values[m][k + 1][k_id-1][t]-cell_values[m][k][k_id-1][t]);
								} else {
									D_1=cell_values[i_TOT - 2 - m][k][k_id-1][t]-cell_values[m][k][k_id-1][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=cell_values[m + 1][k][k_id-1][t]-cell_values[m][k][k_id-1][t];
							}
							if (k != kMax[0] - 2) {
								D_3=cell_values[m][k + 1][k_id-1][t]-cell_values[m][k][k_id-1][t];
							}
							if (m != 0) {
								D_4=cell_values[m - 1][k][k_id-1][t]-cell_values[m][k][k_id-1][t];
							}
						aux_values[m][k][12][t]=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
						aux_values[m][k][13][t]=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						break;
				}
				if (grad_scheme!=1) {	//Vorticity and Lamb vector are calculated in subsequent stages in case the FD gradient scheme is selected
					//Vorticity
					cell_values[m][k][numVars+compVarsN+12][t]=-cell_values[m][k][numVars+compVarsN+4][t]+cell_values[m][k][numVars+compVarsN+3][t];
					//Lamb vector (re-defined in the case of RANS fields)
					cell_values[m][k][numVars+compVarsN+13][t]=cell_values[m][k][vz_id-1][t]*cell_values[m][k][numVars+compVarsN+12][t];
					cell_values[m][k][numVars+compVarsN+14][t]=-cell_values[m][k][vx_id-1][t]*cell_values[m][k][numVars+compVarsN+12][t];
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of computing additional variables at cells centres.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Transfer some variables from cells centres to grid nodes and compute gradients at grid nodes in case of FD scheme selected.
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		for (int k = 0; k <kMax[0]; k++) {
			for (int m = 0; m <i_TOT; m++) {
				//TRANSFER density gradient, Lamb vector, velocity gradient and turbulent kinetic energy gradient
				if (grad_scheme!=1) {		//Only transfer these variables in case they are not directly computed at grid nodes (grad_scheme!=1)
					double inv_d1 = 0;
					double inv_d2 = 0;
					double inv_d3 = 0;
					double inv_d4 = 0;
					if (k==kMax[0]-1) {	//Node is on far-field boundary (C-edge)
						if (m==0) {
							aux_values_N[m][k][17][t] = cell_values[m][k - 1][numVars + compVarsN + 8][t];
							aux_values_N[m][k][18][t] = cell_values[m][k - 1][numVars + compVarsN + 9][t];
							aux_values_N[m][k][19][t] = cell_values[m][k - 1][numVars + compVarsN + 13][t];
							aux_values_N[m][k][20][t] = cell_values[m][k - 1][numVars + compVarsN + 14][t];
							aux_values_N[m][k][21][t] = cell_values[m][k - 1][numVars + compVarsN + 2][t];
							aux_values_N[m][k][22][t] = cell_values[m][k - 1][numVars + compVarsN + 3][t];
							aux_values_N[m][k][23][t] = cell_values[m][k - 1][numVars + compVarsN + 4][t];
							aux_values_N[m][k][24][t] = cell_values[m][k - 1][numVars + compVarsN + 5][t];
							aux_values_N[m][k][30][t] = aux_values[m][k - 1][12][t];
							aux_values_N[m][k][31][t] = aux_values[m][k - 1][13][t];
						} else if (m==i_TOT-1) {
							aux_values_N[m][k][17][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 8][t];
							aux_values_N[m][k][18][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 9][t];
							aux_values_N[m][k][19][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 13][t];
							aux_values_N[m][k][20][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 14][t];
							aux_values_N[m][k][21][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 2][t];
							aux_values_N[m][k][22][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 3][t];
							aux_values_N[m][k][23][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 4][t];
							aux_values_N[m][k][24][t] = cell_values[m - 1][k - 1][numVars + compVarsN + 5][t];
							aux_values_N[m][k][30][t] = aux_values[m - 1][k - 1][12][t];
							aux_values_N[m][k][31][t] = aux_values[m - 1][k - 1][13][t];
						} else {
							inv_d1= 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d2= 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][17][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 8][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 8][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][18][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 9][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 9][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][19][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 13][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 13][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][20][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 14][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 14][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][21][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 2][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 2][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][22][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 3][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 3][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][23][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 4][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 4][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][24][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 5][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 5][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][30][t] = (aux_values[m - 1][k - 1][12][t] * inv_d1 + aux_values[m][k - 1][12][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][31][t] = (aux_values[m - 1][k - 1][13][t] * inv_d1 + aux_values[m][k - 1][13][t] * inv_d2) / (inv_d1 + inv_d2);
						}
					} else if (k==0) {
						if (m==0) {
							inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[i_TOT - 2][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[i_TOT - 2][k][numVars + compVarsN + 16][t]- Z_POLE[t], 2)));
							inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][17][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 8][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 8][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][18][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 9][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 9][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][19][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 13][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 13][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][20][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 14][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 14][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][21][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 2][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 2][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][22][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 3][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 3][t]* inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][23][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 4][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 4][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][24][t] = (cell_values[i_TOT - 2][k][numVars + compVarsN + 5][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 5][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][30][t] = (aux_values[i_TOT - 2][k][12][t] * inv_d2 + aux_values[m][k][12][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][31][t] = (aux_values[i_TOT - 2][k][13][t] * inv_d2 + aux_values[m][k][13][t] * inv_d3) / (inv_d2 + inv_d3);
						} else if (m==i_TOT-1) {
							inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[0][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[0][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][17][t] = (cell_values[m - 1][k][numVars + compVarsN + 8][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 8][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][18][t] = (cell_values[m - 1][k][numVars + compVarsN + 9][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 9][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][19][t] = (cell_values[m - 1][k][numVars + compVarsN + 13][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 13][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][20][t] = (cell_values[m - 1][k][numVars + compVarsN + 14][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 14][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][21][t] = (cell_values[m - 1][k][numVars + compVarsN + 2][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 2][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][22][t] = (cell_values[m - 1][k][numVars + compVarsN + 3][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 3][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][23][t] = (cell_values[m - 1][k][numVars + compVarsN + 4][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 4][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][24][t] = (cell_values[m - 1][k][numVars + compVarsN + 5][t] * inv_d4 + cell_values[0][k][numVars + compVarsN + 5][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][30][t] = (aux_values[m - 1][k][12][t] * inv_d4 + aux_values[0][k][12][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][31][t] = (aux_values[m - 1][k][13][t] * inv_d4 + aux_values[0][k][13][t] * inv_d1) / (inv_d1 + inv_d4);
						} else if (near_ID[m]) {	//Node is on body surface
							inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][17][t] = (cell_values[m - 1][k][numVars + compVarsN + 8][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 8][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][18][t] = (cell_values[m - 1][k][numVars + compVarsN + 9][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 9][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][19][t] = (cell_values[m - 1][k][numVars + compVarsN + 13][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 13][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][20][t] = (cell_values[m - 1][k][numVars + compVarsN + 14][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 14][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][21][t] = (cell_values[m - 1][k][numVars + compVarsN + 2][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 2][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][22][t] = (cell_values[m - 1][k][numVars + compVarsN + 3][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 3][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][23][t] = (cell_values[m - 1][k][numVars + compVarsN + 4][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 4][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][24][t] = (cell_values[m - 1][k][numVars + compVarsN + 5][t] * inv_d4 + cell_values[m][k][numVars + compVarsN + 5][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][30][t] = (aux_values[m - 1][k][12][t] * inv_d4 + aux_values[m][k][12][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][31][t] = (aux_values[m - 1][k][13][t] * inv_d4 + aux_values[m][k][13][t] * inv_d3) / (inv_d3 + inv_d4);
						} else {	//Node is inside fluid domain
							inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][17][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 8][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 8][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 8][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 8][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][18][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 9][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 9][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 9][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 9][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][19][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 13][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 13][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 13][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 13][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][20][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 14][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 14][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 14][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 14][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][21][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 2][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 2][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 2][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 2][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][22][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 3][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 3][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 3][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 3][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][23][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 4][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 4][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 4][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 4][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][24][t] = (cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 5][t] * inv_d1 + cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 5][t] * inv_d2 +
																cell_values[m][k][numVars + compVarsN + 5][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 5][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][30][t] = (aux_values[i_TOT - 1 - m][k][12][t] * inv_d1 + aux_values[i_TOT - 1 - m - 1][k][12][t] * inv_d2 +
																aux_values[m][k][12][t] * inv_d3 + aux_values[m - 1][k][12][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][31][t] = (aux_values[i_TOT - 1 - m][k][13][t] * inv_d1 + aux_values[i_TOT - 1 - m - 1][k][13][t] * inv_d2 +
																aux_values[m][k][13][t] * inv_d3 + aux_values[m - 1][k][13][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
						}
					} else if (m==0) {	//Node is on far-field boundary (outflow edge)
						inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						aux_values_N[m][k][17][t] = (cell_values[m][k - 1][numVars + compVarsN + 8][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 8][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][18][t] = (cell_values[m][k - 1][numVars + compVarsN + 9][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 9][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][19][t] = (cell_values[m][k - 1][numVars + compVarsN + 13][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 13][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][20][t] = (cell_values[m][k - 1][numVars + compVarsN + 14][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 14][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][21][t] = (cell_values[m][k - 1][numVars + compVarsN + 2][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 2][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][22][t] = (cell_values[m][k - 1][numVars + compVarsN + 3][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 3][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][23][t] = (cell_values[m][k - 1][numVars + compVarsN + 4][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 4][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][24][t] = (cell_values[m][k - 1][numVars + compVarsN + 5][t] * inv_d2 + cell_values[m][k][numVars + compVarsN + 5][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][30][t] = (aux_values[m][k - 1][12][t] * inv_d2 + aux_values[m][k][12][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][31][t] = (aux_values[m][k - 1][13][t] * inv_d2 + aux_values[m][k][13][t] * inv_d3) / (inv_d2 + inv_d3);
					} else if (m==i_TOT-1) {
						inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						aux_values_N[m][k][17][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 8][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 8][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][18][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 9][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 9][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][19][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 13][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 13][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][20][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 14][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 14][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][21][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 2][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 2][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][22][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 3][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 3][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][23][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 4][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 4][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][24][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 5][t] * inv_d1 + cell_values[m - 1][k][numVars + compVarsN + 5][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][30][t] = (aux_values[m - 1][k - 1][12][t] * inv_d1 + aux_values[m - 1][k][12][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][31][t] = (aux_values[m - 1][k - 1][13][t] * inv_d1 + aux_values[m - 1][k][13][t] * inv_d4) / (inv_d1 + inv_d4);
					} else {	//Node is in the mesh interior (not on grid boundaries)
					inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
										pow(nodal_values[m][k][2][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
					inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
										pow(nodal_values[m][k][2][t] - cell_values[m][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
					inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
										pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
					inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
										pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
					aux_values_N[m][k][17][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 8][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 8][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 8][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 8][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][18][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 9][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 9][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 9][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 9][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][19][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 13][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 13][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 13][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 13][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][20][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 14][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 14][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 14][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 14][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][21][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 2][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 2][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 2][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 2][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][22][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 3][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 3][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 3][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 3][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][23][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 4][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 4][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 4][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 4][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][24][t] = (cell_values[m - 1][k - 1][numVars + compVarsN + 5][t] * inv_d1 + cell_values[m][k - 1][numVars + compVarsN + 5][t] * inv_d2 +
														cell_values[m][k][numVars + compVarsN + 5][t] * inv_d3 + cell_values[m - 1][k][numVars + compVarsN + 5][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][30][t] = (aux_values[m - 1][k - 1][12][t] * inv_d1 + aux_values[m][k - 1][12][t] * inv_d2 +
														aux_values[m][k][12][t] * inv_d3 + aux_values[m - 1][k][12][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					aux_values_N[m][k][31][t] = (aux_values[m - 1][k - 1][13][t] * inv_d1 + aux_values[m][k - 1][13][t] * inv_d2 +
														aux_values[m][k][13][t] * inv_d3 + aux_values[m - 1][k][13][t] * inv_d4) /
														(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					}
				}
				//COMPUTE gradients of mass density, fluid velocity, 0.5*V^2 and turb. kin. en. (to be implemented), together with Lamb vector, at grid nodes in the case of FD gradient scheme selection
				if (grad_scheme==1) {
					//Local Jacobian
					double J11=0, J12=0, J21=0, J22=0, ddm=0, ddk=0;
					if (m==0) {
						J11=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m+1][k][0][t]-0.5*nodal_values[m+2][k][0][t];
						J12=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m+1][k][2][t]-0.5*nodal_values[m+2][k][2][t];
					} else if (m==i_TOT-1) {
						J11=0.5*nodal_values[m-2][k][0][t]-2.0*nodal_values[m-1][k][0][t]+1.5*nodal_values[m][k][0][t];
						J12=0.5*nodal_values[m-2][k][2][t]-2.0*nodal_values[m-1][k][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J11=0.5*nodal_values[m+1][k][0][t]-0.5*nodal_values[m-1][k][0][t];
                    	J12=0.5*nodal_values[m+1][k][2][t]-0.5*nodal_values[m-1][k][2][t];
					}
                    if (k==0) {
						J21=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k+2][0][t];
						J22=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k+2][2][t];
					} else if (k==kMax[0]-1) {
						J21=0.5*nodal_values[m][k-2][0][t]-2.0*nodal_values[m][k-1][0][t]+1.5*nodal_values[m][k][0][t];
						J22=0.5*nodal_values[m][k-2][2][t]-2.0*nodal_values[m][k-1][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J21=0.5*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k-1][0][t];
						J22=0.5*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k-1][2][t];
					}
                    double detJ=J11*J22-J12*J21;
					//Mass density index gradients
					if (m==0) {
						ddm=-1.5*nodal_values[m][k][rho_id-1][t]+2.0*nodal_values[m+1][k][rho_id-1][t]-0.5*nodal_values[m+2][k][rho_id-1][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*nodal_values[m-2][k][rho_id-1][t]-2.0*nodal_values[m-1][k][rho_id-1][t]+1.5*nodal_values[m][k][rho_id-1][t];
					} else {
						ddm=0.5*nodal_values[m+1][k][rho_id-1][t]-0.5*nodal_values[m-1][k][rho_id-1][t];
					}
					if (k==0) {
						ddk=-1.5*nodal_values[m][k][rho_id-1][t]+2.0*nodal_values[m][k+1][rho_id-1][t]-0.5*nodal_values[m][k+2][rho_id-1][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*nodal_values[m][k-2][rho_id-1][t]-2.0*nodal_values[m][k-1][rho_id-1][t]+1.5*nodal_values[m][k][rho_id-1][t];
					} else {
						ddk=0.5*nodal_values[m][k+1][rho_id-1][t]-0.5*nodal_values[m][k-1][rho_id-1][t];
					}
					//Mass density spatial gradients
					aux_values_N[m][k][17][t] = (J22*ddm - J12*ddk) / detJ;
					aux_values_N[m][k][18][t] = (-J21*ddm + J11*ddk) / detJ;
					//Fluid velocity (x-component) index gradients
					if (m==0) {
						ddm=-1.5*nodal_values[m][k][vx_id-1][t]+2.0*nodal_values[m+1][k][vx_id-1][t]-0.5*nodal_values[m+2][k][vx_id-1][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*nodal_values[m-2][k][vx_id-1][t]-2.0*nodal_values[m-1][k][vx_id-1][t]+1.5*nodal_values[m][k][vx_id-1][t];
					} else {
						ddm=0.5*nodal_values[m+1][k][vx_id-1][t]-0.5*nodal_values[m-1][k][vx_id-1][t];
					}
					if (k==0) {
						ddk=-1.5*nodal_values[m][k][vx_id-1][t]+2.0*nodal_values[m][k+1][vx_id-1][t]-0.5*nodal_values[m][k+2][vx_id-1][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*nodal_values[m][k-2][vx_id-1][t]-2.0*nodal_values[m][k-1][vx_id-1][t]+1.5*nodal_values[m][k][vx_id-1][t];
					} else {
						ddk=0.5*nodal_values[m][k+1][vx_id-1][t]-0.5*nodal_values[m][k-1][vx_id-1][t];
					}
					//Fluid velocity (x-component) spatial gradients
					aux_values_N[m][k][21][t] = (J22*ddm - J12*ddk) / detJ;
					aux_values_N[m][k][22][t] = (-J21*ddm + J11*ddk) / detJ;
					//Fluid velocity (z-component) index gradients
					if (m==0) {
						ddm=-1.5*nodal_values[m][k][vz_id-1][t]+2.0*nodal_values[m+1][k][vz_id-1][t]-0.5*nodal_values[m+2][k][vz_id-1][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*nodal_values[m-2][k][vz_id-1][t]-2.0*nodal_values[m-1][k][vz_id-1][t]+1.5*nodal_values[m][k][vz_id-1][t];
					} else {
						ddm=0.5*nodal_values[m+1][k][vz_id-1][t]-0.5*nodal_values[m-1][k][vz_id-1][t];
					}
					if (k==0) {
						ddk=-1.5*nodal_values[m][k][vz_id-1][t]+2.0*nodal_values[m][k+1][vz_id-1][t]-0.5*nodal_values[m][k+2][vz_id-1][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*nodal_values[m][k-2][vz_id-1][t]-2.0*nodal_values[m][k-1][vz_id-1][t]+1.5*nodal_values[m][k][vz_id-1][t];
					} else {
						ddk=0.5*nodal_values[m][k+1][vz_id-1][t]-0.5*nodal_values[m][k-1][vz_id-1][t];
					}
					//Fluid velocity (z-component) spatial gradients
					aux_values_N[m][k][23][t] = (J22*ddm - J12*ddk) / detJ;
					aux_values_N[m][k][24][t] = (-J21*ddm + J11*ddk) / detJ;
					//Specific kinetic energy per unit mass
					//*************  TEMPORARILY store in aux_values_N[m][k][0][t] and aux_values_N[m][k][1][t] memory locations.
					//*************  They will be overwritten to store other values, once these gradient values will be averaged and stored to cells centres, where they are needed.
					//--- index gradients
					if (m==0) {
						ddm=-1.5*nodal_values[m][k][numVars+1][t]+2.0*nodal_values[m+1][k][numVars+1][t]-0.5*nodal_values[m+2][k][numVars+1][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*nodal_values[m-2][k][numVars+1][t]-2.0*nodal_values[m-1][k][numVars+1][t]+1.5*nodal_values[m][k][numVars+1][t];
					} else {
						ddm=0.5*nodal_values[m+1][k][numVars+1][t]-0.5*nodal_values[m-1][k][numVars+1][t];
					}
					if (k==0) {
						ddk=-1.5*nodal_values[m][k][numVars+1][t]+2.0*nodal_values[m][k+1][numVars+1][t]-0.5*nodal_values[m][k+2][numVars+1][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*nodal_values[m][k-2][numVars+1][t]-2.0*nodal_values[m][k-1][numVars+1][t]+1.5*nodal_values[m][k][numVars+1][t];
					} else {
						ddk=0.5*nodal_values[m][k+1][numVars+1][t]-0.5*nodal_values[m][k-1][numVars+1][t];
					}
					//--- spatial gradients
					aux_values_N[m][k][0][t] = (J22*ddm - J12*ddk) / detJ;
					aux_values_N[m][k][1][t] = (-J21*ddm + J11*ddk) / detJ;
					//Specific turbulent kinetic energy per unit mass
					//--- index gradients
					if (m==0) {
						ddm=-1.5*nodal_values[m][k][k_id-1][t]+2.0*nodal_values[m+1][k][k_id-1][t]-0.5*nodal_values[m+2][k][k_id-1][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*nodal_values[m-2][k][k_id-1][t]-2.0*nodal_values[m-1][k][k_id-1][t]+1.5*nodal_values[m][k][k_id-1][t];
					} else {
						ddm=0.5*nodal_values[m+1][k][k_id-1][t]-0.5*nodal_values[m-1][k][k_id-1][t];
					}
					if (k==0) {
						ddk=-1.5*nodal_values[m][k][k_id-1][t]+2.0*nodal_values[m][k+1][k_id-1][t]-0.5*nodal_values[m][k+2][k_id-1][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*nodal_values[m][k-2][k_id-1][t]-2.0*nodal_values[m][k-1][k_id-1][t]+1.5*nodal_values[m][k][k_id-1][t];
					} else {
						ddk=0.5*nodal_values[m][k+1][k_id-1][t]-0.5*nodal_values[m][k-1][k_id-1][t];
					}
					//--- spatial gradients
					aux_values_N[m][k][30][t] = (J22*ddm - J12*ddk) / detJ;
					aux_values_N[m][k][31][t] = (-J21*ddm + J11*ddk) / detJ;
					//Lamb vector
					aux_values_N[m][k][19][t] =   (aux_values_N[m][k][22][t] - aux_values_N[m][k][23][t]) * nodal_values[m][k][vz_id-1][t];
					aux_values_N[m][k][20][t] = - (aux_values_N[m][k][22][t] - aux_values_N[m][k][23][t]) * nodal_values[m][k][vx_id-1][t];
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of transfering some variables from cells centres to grid nodes.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Average values of nodal gradients and Lamb vector at cells centres, in the case FD gradient method is selected. Vorticicity is computed from averaged velocity gradients.
	start = chrono::steady_clock::now();
	if (grad_scheme==1) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0] - 1; k++) {
				for (int m = 0; m < i_TOT - 1; m++) {
					//Fluid velocity gradient
					cell_values[m][k][numVars+compVarsN+2][t] =
						(aux_values_N[m][k][21][t] + aux_values_N[m + 1][k][21][t] + aux_values_N[m][k + 1][21][t] + aux_values_N[m + 1][k + 1][21][t]) / 4.0;
					cell_values[m][k][numVars+compVarsN+3][t] =
						(aux_values_N[m][k][22][t] + aux_values_N[m + 1][k][22][t] + aux_values_N[m][k + 1][22][t] + aux_values_N[m + 1][k + 1][22][t]) / 4.0;
					cell_values[m][k][numVars+compVarsN+4][t] =
						(aux_values_N[m][k][23][t] + aux_values_N[m + 1][k][23][t] + aux_values_N[m][k + 1][23][t] + aux_values_N[m + 1][k + 1][23][t]) / 4.0;
					cell_values[m][k][numVars+compVarsN+5][t] =
						(aux_values_N[m][k][24][t] + aux_values_N[m + 1][k][24][t] + aux_values_N[m][k + 1][24][t] + aux_values_N[m + 1][k + 1][24][t]) / 4.0;
					//Mass density gradient
					cell_values[m][k][numVars+compVarsN+8][t] =
						(aux_values_N[m][k][17][t] + aux_values_N[m + 1][k][17][t] + aux_values_N[m][k + 1][17][t] + aux_values_N[m + 1][k + 1][17][t]) / 4.0;
					cell_values[m][k][numVars+compVarsN+9][t] =
						(aux_values_N[m][k][18][t] + aux_values_N[m + 1][k][18][t] + aux_values_N[m][k + 1][18][t] + aux_values_N[m + 1][k + 1][18][t]) / 4.0;
					//Specific kinetic energy per unit mass gradient
					cell_values[m][k][numVars+compVarsN+10][t] =
						(aux_values_N[m][k][0][t] + aux_values_N[m + 1][k][0][t] + aux_values_N[m][k + 1][0][t] + aux_values_N[m + 1][k + 1][0][t]) / 4.0;
					cell_values[m][k][numVars+compVarsN+11][t] =
						(aux_values_N[m][k][1][t] + aux_values_N[m + 1][k][1][t] + aux_values_N[m][k + 1][1][t] + aux_values_N[m + 1][k + 1][1][t]) / 4.0;
					//Specific turbulent kinetic energy per unit mass gradient
					aux_values[m][k][12][t] =
						(aux_values_N[m][k][30][t] + aux_values_N[m + 1][k][30][t] + aux_values_N[m][k + 1][30][t] + aux_values_N[m + 1][k + 1][30][t]) / 4.0;
					aux_values[m][k][13][t] =
						(aux_values_N[m][k][31][t] + aux_values_N[m + 1][k][31][t] + aux_values_N[m][k + 1][31][t] + aux_values_N[m + 1][k + 1][31][t]) / 4.0;
					//Vorticity
					cell_values[m][k][numVars+compVarsN+12][t]=-cell_values[m][k][numVars+compVarsN+4][t]+cell_values[m][k][numVars+compVarsN+3][t];
					//Lamb vector
					cell_values[m][k][numVars+compVarsN+13][t] =
						(aux_values_N[m][k][19][t] + aux_values_N[m + 1][k][19][t] + aux_values_N[m][k + 1][19][t] + aux_values_N[m + 1][k + 1][19][t]) / 4.0;
					cell_values[m][k][numVars+compVarsN+14][t] =
						(aux_values_N[m][k][20][t] + aux_values_N[m + 1][k][20][t] + aux_values_N[m][k + 1][20][t] + aux_values_N[m + 1][k + 1][20][t]) / 4.0;
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of averaging values of nodal gradients and Lamb vector at cells centres.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Compute flag values and auxiliar variables at CELLS CENTRES
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		for (int k = 0; k <kMax[0]-1; k++) {
			for (int m = 0; m <i_TOT-1; m++) {
				//Sensors and cells flagging
				//---Shock function (normal Mach number)
				if ((cell_values[m][k][numVars+compVarsN+6][t]!=0)||(cell_values[m][k][numVars+compVarsN+7][t]!=0)) {
					cell_values[m][k][numVars+compVarsN+17][t]=(cell_values[m][k][vx_id-1][t]*cell_values[m][k][numVars+compVarsN+6][t]+
						cell_values[m][k][vz_id-1][t]*cell_values[m][k][numVars+compVarsN+7][t])/(sqrt(gamma_air*R_air*cell_values[m][k][T_id-1][t])*
						sqrt(pow(cell_values[m][k][numVars+compVarsN+6][t],2)+pow(cell_values[m][k][numVars+compVarsN+7][t],2)));
				} else {
					cell_values[m][k][numVars+compVarsN+17][t]=0;
				}
				//---Shock flag (with margins)
				bool check_f0 = false;
				if (flag_id == -1) {
					check_f0 = cell_values[m][k][numVars+compVarsN+17][t]>=Fshock_tsh;
				} else {
					//check_f0 = cell_values[m][k][flag_id-1][t] > 1.5;
					check_f0 = cell_values[m][k][flag_id-1][t] == 2;
				}
				if (check_f0) {
					flags[m][k][0][t]=true;
					for (int kk = k; kk <= k+shock_margin; kk++) {
						if (kk>=0) {
							for (int mm = m-shock_margin; mm <= m+shock_margin; mm++) {
								if ((mm>=0)&&(mm<i_TOT-1)&&(kk<kMax[0]-1)) {
									flags[mm][kk][0][t]=true;
								}
							}
						}
					}
				}
				//------------TEST---------------
				// if ((cell_values[m][k][0][t]>1)&&(cell_values[m][k][0][t]<90)&&(cell_values[m][k][2][t]>=-20)&&(cell_values[m][k][2][t]<=20)) {
				// 	flags[m][k][0][t]=true;
				// }
				//-------------------------------
				if (dom_decomp) {
					//---Eddy viscosity BL sensor
					cell_values[m][k][numVars+compVarsN+18][t]=(cell_values[m][k][numVars+compVarsN+1][t]+cell_values[m][k][mut_id-1][t])/cell_values[m][k][numVars+compVarsN+1][t];
					//---Turbulence dissipation frequency (omega) BL sensor
					cell_values[m][k][numVars+compVarsN+19][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][k_id-1][t]/max(cell_values[m][k][mut_id-1][t],1e-15);
					//---BL flag (with margins)				
					if ((m==0)&&(k==0)) {		//Compute Eddy viscosity BL sensor at freestream only first time cells loop is accessed
						BL_sens_inf = (((0.000001458*sqrt(cell_values[i_TOT/2][kMax[0]-3][T_id-1][t]))/(1+110.4/cell_values[i_TOT/2][kMax[0]-3][T_id-1][t])) + cell_values[i_TOT/2][kMax[0] - 3][mut_id - 1][t]) /
							((0.000001458*sqrt(cell_values[i_TOT/2][kMax[0]-3][T_id-1][t]))/(1+110.4/cell_values[i_TOT/2][kMax[0]-3][T_id-1][t]));
						if (debug_mode) {
							cout << "-----BL_sens_inf at time step " << t << " is: " << BL_sens_inf << "\n";
						}
					}
					bool check_f1 = false;
					if (flag_id == -1) {
						check_f1 = ((cell_values[m][k][numVars+compVarsN+19][t]>=omega_tsh)&&(k!=kMax[0]-2)) || (cell_values[m][k][numVars+compVarsN+18][t]>=K_bl* BL_sens_inf);
					} else {
						//check_f1 = (cell_values[m][k][flag_id-1][t] > 0.5) && (cell_values[m][k][flag_id-1][t] <= 1.5);
						check_f1 = cell_values[m][k][flag_id-1][t] == 1;
					}
					if (check_f1) {
						flags[m][k][1][t]=true;
						int bl_margin_kmax=0;
						int bl_margin_kmin=0;
							if (m<m_te) {
								bl_margin_kmin=bl_margin_NF+((bl_margin_FF_UP-bl_margin_NF)/(cell_values[0][0][0][t]-cell_values[m_te][0][0][t]))*(cell_values[m][k][0][t]-cell_values[m_te][0][0][t]);
								bl_margin_kmax=bl_margin_NF+((bl_margin_FF_LOW-bl_margin_NF)/(cell_values[0][0][0][t]-cell_values[m_te][0][0][t]))*(cell_values[m][k][0][t]-cell_values[m_te][0][0][t]);
							} else {
								if (m>i_TOT - 2 - m_te) {
									bl_margin_kmin=bl_margin_NF+((bl_margin_FF_LOW-bl_margin_NF)/(cell_values[0][0][0][t]-cell_values[m_te][0][0][t]))*(cell_values[m][k][0][t]-cell_values[m_te][0][0][t]);
									bl_margin_kmax=bl_margin_NF+((bl_margin_FF_UP-bl_margin_NF)/(cell_values[0][0][0][t]-cell_values[m_te][0][0][t]))*(cell_values[m][k][0][t]-cell_values[m_te][0][0][t]);
								} else {
									bl_margin_kmin = 0;
									bl_margin_kmax = bl_margin_NF;
								}
							}
						for (int kk = k-bl_margin_kmin; kk <= k+bl_margin_kmax; kk++) {
							if ((kk<kMax[0]-1)&&(-kk-1<kMax[0]-1)) { //&&(cell_values[m][k][numVars+compVarsN+17] < Fshock_tsh)) {
								if (kk>=0) {
									flags[m][kk][1][t]=true;
								} else {
									if ((m < m_te) || (m > i_TOT - 2 - m_te)) {
										flags[i_TOT - 2 - m][-kk-1][1][t]=true;
									}
								}
							}
						}
						if ((m >= m_te)&&(m <= i_TOT - 2 - m_te)) {	//apply wall margin only around body
							for (int kk = 0; kk < wall_margin; kk++) {
								flags[m][kk][1][t]=true;
							}
						}
					}
					//---Non-dimensional entropy variation (w.r.t. freestream) sensor
					cell_values[m][k][numVars+compVarsN+20][t]=log((rho_inf/cell_values[m][k][rho_id-1][t])*pow(cell_values[m][k][T_id-1][t]/T_inf,1/(gamma_air-1)));
					//---Shock wake flag
					if ((cell_values[m][k][numVars+compVarsN+20][t]>=entr_tsh)&&(m>=m_ds)&&(m<=i_TOT-2-m_ds)) {	//&&(!flags[m][k][1])) {
						flags[m][k][2][t]=true;
						flags[m][k][3][t]=true;
					}
					bool check_f0 = false;
					if (flag_id == -1) {
						check_f0 = cell_values[m][k][numVars+compVarsN+17][t]>=Fshock_tsh;
					} else {
						//check_f0 = cell_values[m][k][flag_id-1][t] > 1.5;
						check_f0 = cell_values[m][k][flag_id-1][t] == 2;
					}
					if (check_f0) {
					//if (flags[m][k][0]&&(!flags[m][k][1])) {
						flags[m][k][2][t]=true;
						flags[m][k][3][t]=true;
						if (!flags[m][k][1][t]) {
							for (int kk = k; kk <= k+shock_margin_drag_kmax; kk++) {
								for (int mm = m-shock_margin_drag_mmin; mm <= m + shock_margin_drag_mmax; mm++) {
									flags[mm][kk][2][t]=true;
									flags[mm][kk][3][t]=true;
								}
							}
						}
					}
				}
				//Auxiliar variables (for easier elementary aerodynamic force calculation: FS_omega terms)
				aux_values[m][k][0][t]=(cell_values[m][k][numVars+compVarsN+15][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t]+
					cell_values[m][k][numVars+compVarsN+16][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t]);
				aux_values[m][k][1][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t];
				aux_values[m][k][2][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t];
				//Viscous stress tensor at cells centres
				if (NS) {
					if (Lamb_form==0) {	//Include Reynolds stresses
						aux_values[m][k][3][t]=2*(cell_values[m][k][numVars+compVarsN+1][t]+cell_values[m][k][mut_id-1][t])*(cell_values[m][k][numVars+compVarsN+2][t]-
											(1.0/3.0)*(cell_values[m][k][numVars+compVarsN+2][t]+cell_values[m][k][numVars+compVarsN+5][t]));		//2*(mu+mu_T)*(gradV)_0S(1,1)
						aux_values[m][k][4][t]=2*(cell_values[m][k][numVars+compVarsN+1][t]+cell_values[m][k][mut_id-1][t])*(cell_values[m][k][numVars+compVarsN+3][t]+
											cell_values[m][k][numVars+compVarsN+4][t])/2;													//2*(mu+mu_T)*(gradV)_0S(1,3)
						aux_values[m][k][5][t]=2*(cell_values[m][k][numVars+compVarsN+1][t]+cell_values[m][k][mut_id-1][t])*(
											-(1.0/3.0)*(cell_values[m][k][numVars+compVarsN+2][t]+cell_values[m][k][numVars+compVarsN+5][t]));	//2*(mu+mu_T)*(gradV)_0S(2,2)
						aux_values[m][k][6][t]=2*(cell_values[m][k][numVars+compVarsN+1][t]+cell_values[m][k][mut_id-1][t])*(cell_values[m][k][numVars+compVarsN+5][t]-
											(1.0/3.0)*(cell_values[m][k][numVars+compVarsN+2][t]+cell_values[m][k][numVars+compVarsN+5][t]));		//2*(mu+mu_T)*(gradV)_0S(3,3)
					} else {	//Don't include Reynolds stresses
						aux_values[m][k][3][t]=2*cell_values[m][k][numVars+compVarsN+1][t]*(cell_values[m][k][numVars+compVarsN+2][t]-
											(1.0/3.0)*(cell_values[m][k][numVars+compVarsN+2][t]+cell_values[m][k][numVars+compVarsN+5][t]));		//2*(mu+mu_T)*(gradV)_0S(1,1)
						aux_values[m][k][4][t]=2*cell_values[m][k][numVars+compVarsN+1][t]*(cell_values[m][k][numVars+compVarsN+3][t]+
											cell_values[m][k][numVars+compVarsN+4][t])/2;													//2*(mu+mu_T)*(gradV)_0S(1,3)
						aux_values[m][k][5][t]=2*cell_values[m][k][numVars+compVarsN+1][t]*(
											-(1.0/3.0)*(cell_values[m][k][numVars+compVarsN+2][t]+cell_values[m][k][numVars+compVarsN+5][t]));	//2*(mu+mu_T)*(gradV)_0S(2,2)
						aux_values[m][k][6][t]=2*cell_values[m][k][numVars+compVarsN+1][t]*(cell_values[m][k][numVars+compVarsN+5][t]-
											(1.0/3.0)*(cell_values[m][k][numVars+compVarsN+2][t]+cell_values[m][k][numVars+compVarsN+5][t]));		//2*(mu+mu_T)*(gradV)_0S(3,3)
					}
				}
				//Auxiliar variables (for easier elementary aerodynamic force calculation: F_grad(rho)_omega terms)
				aux_values[m][k][9][t]=(cell_values[m][k][numVars+compVarsN+15][t]*cell_values[m][k][numVars+1][t]*cell_values[m][k][numVars+compVarsN+8][t]+
					cell_values[m][k][numVars+compVarsN+16][t]*cell_values[m][k][numVars+1][t]*cell_values[m][k][numVars+compVarsN+9][t]);
				aux_values[m][k][10][t]=cell_values[m][k][numVars+1][t]*cell_values[m][k][numVars+compVarsN+8][t];
				aux_values[m][k][11][t]=cell_values[m][k][numVars+1][t]*cell_values[m][k][numVars+compVarsN+9][t];
			}
		}
		//---Spatial loops closed and re-opened below because shock flag status of grid cells might have been modified when flagging other cells due to margins application
		if (dom_decomp) {
			for (int k = 0; k <kMax[0]-1; k++) {		//Separation of viscous and wave regions
				for (int m = 0; m <i_TOT-1; m++) {
					// if (flags[m][k][0]) {
					// 	flags[m][k][2] = true;
					// 	flags[m][k][3] = true;
					// }
					if (NS) {	//Serapate viscous and wave drag regions, with treshold at shock foot
						bool check_f = false;
						if (flag_id == -1) {
							check_f = cell_values[m][k][numVars+compVarsN+17][t] >= Fshock_tail;
						} else {
							//check_f = cell_values[m][k][flag_id-1][t] == 2;
							check_f = cell_values[m][k][numVars+compVarsN+17][t] >= Fshock_tail;
						}
						if (flags[m][k][1][t] && (!check_f)) {
							flags[m][k][2][t] = false;
							flags[m][k][3][t] = false;
						} else if (flags[m][k][1][t] && (check_f)) {
							flags[m][k][1][t] = false;
						}
					} else {	//Set NO VISCOUS REGION in case of effectively inviscid flow computations
						flags[m][k][1][t] = false;
					}
				}
			}
			if (debug_mode) {
				for (int k = 0; k <kMax[0]-1; k++) {
					for (int m = 0; m <i_TOT-1; m++) {
						flag_org[m][k][t]=flags[m][k][2][t];
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of computing flags and auxiliary variables at cells centres.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Compute auxiliary nodal variables (for easier elementary aerodynamic force calculation)
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		for (int k = 0; k <kMax[0]; k++) {
			for (int m = 0; m <i_TOT; m++) {
				//COMPUTE auxiliary nodal variables (FS_Omega terms)
				if (!flux_scheme) {
					aux_values_N[m][k][0][t]=((nodal_values[m][k][0][t] - X_POLE[t]) * nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][19][t] +
											(nodal_values[m][k][2][t] - Z_POLE[t]) * nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][20][t]);
					aux_values_N[m][k][1][t]=nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][19][t];
					aux_values_N[m][k][2][t]=nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][20][t];
				}
				//COMPUTE variables used in thermodynamic methods
				if (thd_methods || (Lamb_form>0) ) {
					//---Non-dimensional entropy variation (w.r.t. freestream)
					if (INC) {
						aux_values_N[m][k][3][t] = log((p_inf / nodal_values[m][k][numVars][t]) * pow(nodal_values[m][k][T_id - 1][t] / T_inf, gamma_air / (gamma_air - 1)));
					} else {
						aux_values_N[m][k][3][t] = log((rho_inf / nodal_values[m][k][rho_id - 1][t]) * pow(nodal_values[m][k][T_id - 1][t] / T_inf, 1 / (gamma_air - 1)));
					}
					//---Non-dimensional total enthalpy variation (w.r.t. freestream)
					aux_values_N[m][k][4][t] = (((gamma_air * R_air / (gamma_air - 1)) * nodal_values[m][k][T_id - 1][t] + (pow(nodal_values[m][k][vx_id - 1][t], 2) + pow(nodal_values[m][k][vz_id - 1][t], 2)) / 2) -
						((gamma_air * R_air / (gamma_air - 1)) * T_inf + pow(V_inf, 2) / 2)) / pow(V_inf, 2);
					//---Non-dimensional pressure variation (w.r.t. freestream)
					aux_values_N[m][k][5][t] = (nodal_values[m][k][numVars][t] - p_inf) / p_inf;
					//---Modified Negative g-function (thermodynamic method - Paparone/Tognaccini)
					aux_values_N[m][k][6][t] = (nodal_values[m][k][vx_id - 1][t] / sqrt(pow(nodal_values[m][k][vx_id - 1][t], 2) + pow(nodal_values[m][k][vz_id - 1][t], 2)))*(
						(-1 / (gamma_air * pow(M_inf, 2))) * aux_values_N[m][k][3][t] -
						((1 + (gamma_air - 1) * pow(M_inf, 2)) / (2 * pow(gamma_air, 2) * pow(M_inf, 4))) * pow(aux_values_N[m][k][3][t], 2));
					if (std::isnan(aux_values_N[m][k][6][t])) {
						aux_values_N[m][k][6][t] = 0;
					}
					//---Modified f-function (thermodynamic method - Paparone/Tognaccini, complete formula)
					//aux_values_N[m][k][7] = (nodal_values[m][k][vx_id-1]/(sqrt(pow(nodal_values[m][k][vx_id - 1],2)+pow(nodal_values[m][k][vz_id - 1],2)))) *
						//sqrt(1 + 2 * aux_values_N[m][k][4] - (2 / ((gamma_air - 1) * pow(M_inf, 2))) *
						//(pow(aux_values_N[m][k][5] + 1, (gamma_air - 1) / gamma_air) * exp(aux_values_N[m][k][3] * (gamma_air - 1) / gamma_air) - 1));
					//---f-function (thermodynamic method - Destarac & VdV)
					aux_values_N[m][k][7][t] =
						sqrt(1 + 2 * aux_values_N[m][k][4][t] - (2 / ((gamma_air - 1) * pow(M_inf, 2))) *
						(exp(aux_values_N[m][k][3][t] * (gamma_air - 1) / gamma_air) - 1));
					if (std::isnan(aux_values_N[m][k][7][t])) {
						aux_values_N[m][k][7][t] = 0;
					}
					//---Pressure variation term (thermodynamic method - Paparone/Tognaccini, complete formula)
					aux_values_N[m][k][8][t] = rho_inf*V_inf * aux_values_N[m][k][5][t] * (1 / (gamma_air * pow(M_inf, 2)));
				}
				//COMPUTE viscous stress tensor at grid nodes
				if (NS) {
					if (Lamb_form==0) {	//Include Reynolds stresses
						aux_values_N[m][k][9][t]=2*(((0.000001458*sqrt(nodal_values[m][k][T_id-1][t]))/(1+110.4/nodal_values[m][k][T_id-1][t]))+nodal_values[m][k][mut_id-1][t])*
								(aux_values_N[m][k][21][t]-(1.0/3.0)*(aux_values_N[m][k][21][t]+aux_values_N[m][k][24][t]));		//2*(mu+mu_T)*(gradV)_0S(1,1)
						aux_values_N[m][k][10][t]=2*(((0.000001458*sqrt(nodal_values[m][k][T_id-1][t]))/(1+110.4/nodal_values[m][k][T_id-1][t]))+nodal_values[m][k][mut_id-1][t])*
								(aux_values_N[m][k][22][t]+aux_values_N[m][k][23][t])/2;													//2*(mu+mu_T)*(gradV)_0S(1,3)
						aux_values_N[m][k][11][t]=2*(((0.000001458*sqrt(nodal_values[m][k][T_id-1][t]))/(1+110.4/nodal_values[m][k][T_id-1][t]))+nodal_values[m][k][mut_id-1][t])*
								(aux_values_N[m][k][24][t]-(1.0/3.0)*(aux_values_N[m][k][21][t]+aux_values_N[m][k][24][t]));		//2*(mu+mu_T)*(gradV)_0S(3,3)
					} else {	//Don't include Reynolds stresses
						aux_values_N[m][k][9][t]=2*((0.000001458*sqrt(nodal_values[m][k][T_id-1][t]))/(1+110.4/nodal_values[m][k][T_id-1][t]))*
								(aux_values_N[m][k][21][t]-(1.0/3.0)*(aux_values_N[m][k][21][t]+aux_values_N[m][k][24][t]));		//2*(mu+mu_T)*(gradV)_0S(1,1)
						aux_values_N[m][k][10][t]=2*((0.000001458*sqrt(nodal_values[m][k][T_id-1][t]))/(1+110.4/nodal_values[m][k][T_id-1][t]))*
								(aux_values_N[m][k][22][t]+aux_values_N[m][k][23][t])/2;													//2*(mu+mu_T)*(gradV)_0S(1,3)
						aux_values_N[m][k][11][t]=2*((0.000001458*sqrt(nodal_values[m][k][T_id-1][t]))/(1+110.4/nodal_values[m][k][T_id-1][t]))*
								(aux_values_N[m][k][24][t]-(1.0/3.0)*(aux_values_N[m][k][21][t]+aux_values_N[m][k][24][t]));		//2*(mu+mu_T)*(gradV)_0S(3,3)
					}
				}
				//COMPUTE auxiliary nodal variables (F_grad(rho) terms)
				aux_values_N[m][k][14][t]=((nodal_values[m][k][0][t] - X_POLE[t]) * nodal_values[m][k][numVars+1][t] * aux_values_N[m][k][17][t] +
										(nodal_values[m][k][2][t] - Z_POLE[t]) * nodal_values[m][k][numVars+1][t] * aux_values_N[m][k][18][t]);
				aux_values_N[m][k][15][t]=nodal_values[m][k][numVars+1][t]*aux_values_N[m][k][17][t];
				aux_values_N[m][k][16][t]=nodal_values[m][k][numVars+1][t]*aux_values_N[m][k][18][t];
				//COMPUTE time derivatives of fluid velocity and mass density
				if (Nt>=3) {
					//---Non-inertial time derivatives of absolute (U_FORM=0) or relative (U_FORM=1) fluid velocity
					if (t==0) {				//2nd order forward time derivative (uniform time steps)
						aux_values_N[m][k][25][t] = -1.5*nodal_values[m][k][vx_id - 1][t] + 2*nodal_values[m][k][vx_id - 1][t+1] - 0.5*nodal_values[m][k][vx_id - 1][t+2];
						aux_values_N[m][k][26][t] = -1.5*nodal_values[m][k][vz_id - 1][t] + 2*nodal_values[m][k][vz_id - 1][t+1] - 0.5*nodal_values[m][k][vz_id - 1][t+2];
					} else if (t==Nt-1) {	//2nd order backward time derivative (uniform time steps)
						aux_values_N[m][k][25][t] = 0.5*nodal_values[m][k][vx_id - 1][t-2] - 2*nodal_values[m][k][vx_id - 1][t-1] + 1.5*nodal_values[m][k][vx_id - 1][t];
						aux_values_N[m][k][26][t] = 0.5*nodal_values[m][k][vz_id - 1][t-2] - 2*nodal_values[m][k][vz_id - 1][t-1] + 1.5*nodal_values[m][k][vz_id - 1][t];
					} else {				//2nd order central time derivative (uniform time steps)
						aux_values_N[m][k][25][t] = - 0.5*nodal_values[m][k][vx_id - 1][t-1] + 0.5*nodal_values[m][k][vx_id - 1][t+1];
						aux_values_N[m][k][26][t] = - 0.5*nodal_values[m][k][vz_id - 1][t-1] + 0.5*nodal_values[m][k][vz_id - 1][t+1];
					}
					aux_values_N[m][k][25][t]=aux_values_N[m][k][25][t]/(time[1]-time[0]);
					aux_values_N[m][k][26][t]=aux_values_N[m][k][26][t]/(time[1]-time[0]);
					if (U_FORM==0) {			//Inertial time derivatives of absolute fluid velocity
						aux_values_N[m][k][27][t]=aux_values_N[m][k][25][t]-
												(nodal_values[m][k][vgx_id-1][t]*aux_values_N[m][k][21][t]+nodal_values[m][k][vgz_id-1][t]*aux_values_N[m][k][22][t])+
												freq_a*ampl_a*cos(freq_a*(time[t]-t0_a))*nodal_values[m][k][vz_id-1][t];
						aux_values_N[m][k][28][t]=aux_values_N[m][k][26][t]-
												(nodal_values[m][k][vgx_id-1][t]*aux_values_N[m][k][23][t]+nodal_values[m][k][vgz_id-1][t]*aux_values_N[m][k][24][t])-
												freq_a*ampl_a*cos(freq_a*(time[t]-t0_a))*nodal_values[m][k][vx_id-1][t];
					} else if (U_FORM==1) {		//Mass specific apparent force of relative motion
						aux_values_N[m][k][27][t] = freq_x*freq_x*ampl_x*sin(freq_x*(time[t]-t0_x)) + 
													freq_a*freq_a*ampl_a*sin(freq_a*(time[t]-t0_a)) * (nodal_values[m][k][2][t]-(z_CoR + ampl_z*sin(freq_z*(time[t]-t0_z)))) +
													freq_a*freq_a*ampl_a*ampl_a*cos(freq_a*(time[t]-t0_a))*cos(freq_a*(time[t]-t0_a)) * (nodal_values[m][k][0][t]-(x_CoR + ampl_x*sin(freq_x*(time[t]-t0_x)))) -
													2*freq_a*ampl_a*cos(freq_a*(time[t]-t0_a)) * nodal_values[m][k][vz_id-1][t];
						aux_values_N[m][k][28][t] = freq_z*freq_z*ampl_z*sin(freq_z*(time[t]-t0_z)) - 
													freq_a*freq_a*ampl_a*sin(freq_a*(time[t]-t0_a)) * (nodal_values[m][k][0][t]-(x_CoR + ampl_x*sin(freq_x*(time[t]-t0_x)))) +
													freq_a*freq_a*ampl_a*ampl_a*cos(freq_a*(time[t]-t0_a))*cos(freq_a*(time[t]-t0_a)) * (nodal_values[m][k][2][t]-(z_CoR + ampl_z*sin(freq_z*(time[t]-t0_z)))) +
													2*freq_a*ampl_a*cos(freq_a*(time[t]-t0_a)) * nodal_values[m][k][vx_id-1][t];
					} else if (U_FORM==2) {		//Non-inertial time derivatives of relative fluid velocity
						if (t==0) {				//2nd order forward time derivative (uniform time steps)
							aux_values_N[m][k][27][t] = -1.5*(nodal_values[m][k][vx_id - 1][t]-nodal_values[m][k][vgx_id - 1][t])+
																2*(nodal_values[m][k][vx_id - 1][t+1]-nodal_values[m][k][vgx_id - 1][t+1]) -
																0.5*(nodal_values[m][k][vx_id - 1][t+2]-nodal_values[m][k][vgx_id - 1][t+2]);
							aux_values_N[m][k][28][t] = -1.5*(nodal_values[m][k][vz_id - 1][t]-nodal_values[m][k][vgz_id - 1][t])+
																2*(nodal_values[m][k][vz_id - 1][t+1]-nodal_values[m][k][vgz_id - 1][t+1]) -
																0.5*(nodal_values[m][k][vz_id - 1][t+2]-nodal_values[m][k][vgz_id - 1][t+2]);
						} else if (t==Nt-1) {	//2nd order backward time derivative (uniform time steps)
							aux_values_N[m][k][27][t] = 0.5*(nodal_values[m][k][vx_id - 1][t-2]-nodal_values[m][k][vgx_id - 1][t-2]) -
																2*(nodal_values[m][k][vx_id - 1][t-1]-nodal_values[m][k][vgx_id - 1][t-1]) +
																1.5*(nodal_values[m][k][vx_id - 1][t]-nodal_values[m][k][vgx_id - 1][t]);
							aux_values_N[m][k][28][t] = 0.5*(nodal_values[m][k][vz_id - 1][t-2]-nodal_values[m][k][vgz_id - 1][t-2]) -
																2*(nodal_values[m][k][vz_id - 1][t-1]-nodal_values[m][k][vgz_id - 1][t-1]) +
																1.5*(nodal_values[m][k][vz_id - 1][t]-nodal_values[m][k][vgz_id - 1][t]);
						} else {				//2nd order central time derivative (uniform time steps)
							aux_values_N[m][k][27][t] = - 0.5*(nodal_values[m][k][vx_id - 1][t-1]-nodal_values[m][k][vgx_id - 1][t-1]) +
																	0.5*(nodal_values[m][k][vx_id - 1][t+1]-nodal_values[m][k][vgx_id - 1][t+1]);
							aux_values_N[m][k][28][t] = - 0.5*(nodal_values[m][k][vz_id - 1][t-1]-nodal_values[m][k][vgz_id - 1][t-1]) +
																	0.5*(nodal_values[m][k][vz_id - 1][t+1]-nodal_values[m][k][vgz_id - 1][t+1]);
						}
						aux_values_N[m][k][27][t]=aux_values_N[m][k][27][t]/(time[1]-time[0]);
						aux_values_N[m][k][28][t]=aux_values_N[m][k][28][t]/(time[1]-time[0]);
					}
					//---Non-inertial time derivatives of mass density
					if (t==0) {				//2nd order forward time derivative (uniform time steps)
						aux_values_N[m][k][29][t] = -1.5*nodal_values[m][k][rho_id - 1][t] + 2*nodal_values[m][k][rho_id - 1][t+1] - 0.5*nodal_values[m][k][rho_id - 1][t+2];
					} else if (t==Nt-1) {	//2nd order backward time derivative (uniform time steps)
						aux_values_N[m][k][29][t] = 0.5*nodal_values[m][k][rho_id - 1][t-2] - 2*nodal_values[m][k][rho_id - 1][t-1] + 1.5*nodal_values[m][k][rho_id - 1][t];
					} else {				//2nd order central time derivative (uniform time steps)
						aux_values_N[m][k][29][t] = - 0.5*nodal_values[m][k][rho_id - 1][t-1] + 0.5*nodal_values[m][k][rho_id - 1][t+1];
					}
					aux_values_N[m][k][29][t]=aux_values_N[m][k][29][t]/(time[1]-time[0]);
					if (U_FORM==0) {	//Convert to inertial time derivatives of mass density
						aux_values_N[m][k][29][t]=aux_values_N[m][k][29][t]-
												(nodal_values[m][k][vgx_id-1][t]*aux_values_N[m][k][17][t]+nodal_values[m][k][vgz_id-1][t]*aux_values_N[m][k][18][t]);
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of computing auxiliary nodal variables.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Average values of some auxiliary nodal variables at cells centres, in the case WLS gradient method is selected and LAMB_FORM > 0.
	//---TEMPORARY auxiliary variables are used (and reused for other quantities afterwards in the code...)
	start = chrono::steady_clock::now();
	if ((grad_scheme==2)&&(Lamb_form>0)) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0] - 1; k++) {
				for (int m = 0; m < i_TOT - 1; m++) {
					//Specific total enthalpy per unit mass (dimensional)
					aux_values[m][k][14][t] = pow(V_inf, 2) *
						(aux_values_N[m][k][4][t] + aux_values_N[m + 1][k][4][t] + aux_values_N[m][k + 1][4][t] + aux_values_N[m + 1][k + 1][4][t]) / 4.0;
					//Specific entropy variation per unit mass  (dimensional)
					aux_values[m][k][15][t] = R_air *
						(aux_values_N[m][k][3][t] + aux_values_N[m + 1][k][3][t] + aux_values_N[m][k + 1][3][t] + aux_values_N[m + 1][k + 1][3][t]) / 4.0;
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of averaging values of some auxiliary nodal variables.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Divergence of the viscous stress tensor
	start = chrono::steady_clock::now();
	if (NS) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k <kMax[0]-1; k++) {
				for (int m = 0; m <i_TOT-1; m++) {
					cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
								rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
					switch (grad_scheme) {
                    	case 0:     //Green-Gauss (nodes-based)
							if (flux_scheme) {
								if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
									aux_values[m][k][7][t] = 0;
								} else if (k != 0) {
									aux_values[m][k][7][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
										dz_1 * (a_1 * aux_values[m][k - 1][3][t] + (1 - a_1) * aux_values[m][k][3][t]) - dx_1 * (a_1 * aux_values[m][k - 1][4][t] + (1 - a_1) * aux_values[m][k][4][t]) +
										dz_2 * (a_2 * aux_values[m + 1][k][3][t] + (1 - a_2) * aux_values[m][k][3][t]) - dx_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) +
										dz_3 * (a_3 * aux_values[m][k + 1][3][t] + (1 - a_3) * aux_values[m][k][3][t]) - dx_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) +
										dz_4 * (a_4 * aux_values[m - 1][k][3][t] + (1 - a_4) * aux_values[m][k][3][t]) - dx_4 * (a_4 * aux_values[m - 1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]));
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface
										if (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option
											aux_values[m][k][7][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
												dz_1 * (a_1 * aux_values[m][k][3][t] + (1 - a_1) * aux_values[m][k][3][t]) - dx_1  * (a_1 * aux_values[m][k][4][t] + (1 - a_1) * aux_values[m][k][4][t]) +
												dz_2 * (a_2 * aux_values[m + 1][k][3][t] + (1 - a_2) * aux_values[m][k][3][t]) - dx_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) +
												dz_3 * (a_3 * aux_values[m][k + 1][3][t] + (1 - a_3) * aux_values[m][k][3][t]) - dx_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) +
												dz_4 * (a_4 * aux_values[m - 1][k][3][t] + (1 - a_4) * aux_values[m][k][3][t]) - dx_4 * (a_4 * aux_values[m - 1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]));
										} else {	//In this case the cell_metrics function has set a_1 to allow extrapolation based on the first two cells.
											aux_values[m][k][7][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
												dz_1  * (-a_1 * aux_values[m][k + 1][3][t] + (1 + a_1) * aux_values[m][k][3][t]) - dx_1  * (-a_1 * aux_values[m][k + 1][4][t] + (1 + a_1) * aux_values[m][k][4][t]) +
												dz_2 * (a_2 * aux_values[m + 1][k][3][t] + (1 - a_2) * aux_values[m][k][3][t]) - dx_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) +
												dz_3 * (a_3 * aux_values[m][k + 1][3][t] + (1 - a_3) * aux_values[m][k][3][t]) - dx_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) +
												dz_4 * (a_4 * aux_values[m - 1][k][3][t] + (1 - a_4) * aux_values[m][k][3][t]) - dx_4 * (a_4 * aux_values[m - 1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]));
										}
									} else {
										aux_values[m][k][7][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
											dz_1 * (a_1 * aux_values[i_TOT - 2 - m][k][3][t] + (1 - a_1) * aux_values[m][k][3][t]) - dx_1 * (a_1 * aux_values[i_TOT - 2 - m][k][4][t] + (1 - a_1) * aux_values[m][k][4][t]) +
											dz_2 * (a_2 * aux_values[m + 1][k][3][t] + (1 - a_2) * aux_values[m][k][3][t]) - dx_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) +
											dz_3 * (a_3 * aux_values[m][k + 1][3][t] + (1 - a_3) * aux_values[m][k][3][t]) - dx_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) +
											dz_4 * (a_4 * aux_values[m - 1][k][3][t] + (1 - a_4) * aux_values[m][k][3][t]) - dx_4 * (a_4 * aux_values[m - 1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]));
									}
								}
								if ((m == 0) || (m == i_TOT - 2) || (k == kMax[0] - 2)) {     // cell is at far-field boundary
									aux_values[m][k][8][t] = 0;
								} else if (k != 0) {
									aux_values[m][k][8][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
										dz_1 * (a_1 * aux_values[m][k - 1][4][t] + (1 - a_1) * aux_values[m][k][4][t]) - dx_1 * (a_1 * aux_values[m][k - 1][6][t] + (1 - a_1) * aux_values[m][k][6][t]) +
										dz_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) - dx_2 * (a_2 * aux_values[m + 1][k][6][t] + (1 - a_2) * aux_values[m][k][6][t]) +
										dz_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) - dx_3 * (a_3 * aux_values[m][k + 1][6][t] + (1 - a_3) * aux_values[m][k][6][t]) +
										dz_4 * (a_4 * aux_values[m - 1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]) - dx_4 * (a_4 * aux_values[m - 1][k][6][t] + (1 - a_4) * aux_values[m][k][6][t]));
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {     //cell is adjacent to the body surface	
										if (!MOV_GRD || ((U_FORM==1)||(U_FORM==2))) {	//In this case the cell_metrics function has set a_1=0.5, so 1st cell extrapolation is the option					
											aux_values[m][k][8][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
												dz_1 * (a_1 * aux_values[m][k][4][t] + (1 - a_1) * aux_values[m][k][4][t]) - dx_1  * (a_1 * aux_values[m][k][6][t] + (1 - a_1) * aux_values[m][k][6][t]) +
												dz_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) - dx_2 * (a_2 * aux_values[m + 1][k][6][t] + (1 - a_2) * aux_values[m][k][6][t]) +
												dz_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) - dx_3 * (a_3 * aux_values[m][k + 1][6][t] + (1 - a_3) * aux_values[m][k][6][t]) +
												dz_4 * (a_4 * aux_values[m - 1][k][4][t]+ (1 - a_4) * aux_values[m][k][4][t]) - dx_4 * (a_4 * aux_values[m - 1][k][6][t] + (1 - a_4) * aux_values[m][k][6][t]));
										} else {	//In this case the cell_metrics function has set a_1 to allow extrapolation based on the first two cells.
											aux_values[m][k][8][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
												dz_1  * (-a_1 * aux_values[m][k + 1][4][t] + (1 + a_1) * aux_values[m][k][4][t]) - dx_1  * (-a_1 * aux_values[m][k + 1][6][t] + (1 + a_1) * aux_values[m][k][6][t]) +
												dz_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) - dx_2 * (a_2 * aux_values[m + 1][k][6][t] + (1 - a_2) * aux_values[m][k][6][t]) +
												dz_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) - dx_3 * (a_3 * aux_values[m][k + 1][6][t] + (1 - a_3) * aux_values[m][k][6][t]) +
												dz_4 * (a_4 * aux_values[m - 1][k][4][t]+ (1 - a_4) * aux_values[m][k][4][t]) - dx_4 * (a_4 * aux_values[m - 1][k][6][t] + (1 - a_4) * aux_values[m][k][6][t]));
										}
									} else {
										aux_values[m][k][8][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
											dz_1 * (a_1 * aux_values[i_TOT - 2 - m][k][4][t] + (1 - a_1) * aux_values[m][k][4][t]) - dx_1 * (a_1 * aux_values[i_TOT - 2 - m][k][6][t] + (1 - a_1) * aux_values[m][k][6][t]) +
											dz_2 * (a_2 * aux_values[m + 1][k][4][t] + (1 - a_2) * aux_values[m][k][4][t]) - dx_2 * (a_2 * aux_values[m + 1][k][6][t] + (1 - a_2) * aux_values[m][k][6][t]) +
											dz_3 * (a_3 * aux_values[m][k + 1][4][t] + (1 - a_3) * aux_values[m][k][4][t]) - dx_3 * (a_3 * aux_values[m][k + 1][6][t] + (1 - a_3) * aux_values[m][k][6][t]) +
											dz_4 * (a_4 * aux_values[m - 1][k][4][t] + (1 - a_4) * aux_values[m][k][4][t]) - dx_4 * (a_4 * aux_values[m - 1][k][6][t] + (1 - a_4) * aux_values[m][k][6][t]));
									}
								}
							} else {				
								aux_values[m][k][7][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
									dz_1 * ((aux_values_N[m][k][9][t]+aux_values_N[m+1][k][9][t])/2) - dx_1 * ((aux_values_N[m][k][10][t]+aux_values_N[m+1][k][10][t])/2) +
									dz_2 * ((aux_values_N[m+1][k][9][t]+aux_values_N[m+1][k+1][9][t])/2) - dx_2 * ((aux_values_N[m+1][k][10][t]+aux_values_N[m+1][k+1][10][t])/2) +
									dz_3 * ((aux_values_N[m+1][k+1][9][t]+aux_values_N[m][k+1][9][t])/2) - dx_3 * ((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2) +
									dz_4 * ((aux_values_N[m][k+1][9][t]+aux_values_N[m][k][9][t])/2) - dx_4 * ((aux_values_N[m][k+1][10][t]+aux_values_N[m][k][10][t])/2));
								aux_values[m][k][8][t] = (1/cell_values[m][k][numVars+compVarsN][t])*(
									dz_1 * ((aux_values_N[m][k][10][t]+aux_values_N[m+1][k][10][t])/2) - dx_1 * ((aux_values_N[m][k][11][t]+aux_values_N[m+1][k][11][t])/2) +
									dz_2 * ((aux_values_N[m+1][k][10][t]+aux_values_N[m+1][k+1][10][t])/2) - dx_2 * ((aux_values_N[m+1][k][11][t]+aux_values_N[m+1][k+1][11][t])/2) +
									dz_3 * ((aux_values_N[m+1][k+1][10][t]+aux_values_N[m][k+1][10][t])/2) - dx_3 * ((aux_values_N[m+1][k+1][11][t]+aux_values_N[m][k+1][11][t])/2) +
									dz_4 * ((aux_values_N[m][k+1][10][t]+aux_values_N[m][k][10][t])/2) - dx_4 * ((aux_values_N[m][k+1][11][t]+aux_values_N[m][k][11][t])/2));
							}
							break;
						case 1:     //FD method (structured grids)
							//Viscous stress tensor divergence computed at grid nodes in the following sections, then averaged to cells centres.
							break;
						case 2:     //Weighted Least Squares method
							double D_1, D_2, D_3, D_4;
							//---tau_xz
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=aux_values[m][k - 1][4][t]-aux_values[m][k][4][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									//no action, D_1 equal zero (initialization) is ok
								} else {
									D_1=aux_values[i_TOT - 2 - m][k][4][t]-aux_values[m][k][4][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=aux_values[m + 1][k][4][t]-aux_values[m][k][4][t];
							}
							if (k != kMax[0] - 2) {
								D_3=aux_values[m][k + 1][4][t]-aux_values[m][k][4][t];
							}
							if (m != 0) {
								D_4=aux_values[m - 1][k][4][t]-aux_values[m][k][4][t];
							}
							aux_values[m][k][8][t] = aux_values[m][k][8][t] + WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
							aux_values[m][k][7][t] = aux_values[m][k][7][t] + WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
							//---tau_xx
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=aux_values[m][k - 1][3][t]-aux_values[m][k][3][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									//no action, D_1 equal zero (initialization) is ok
								} else {
									D_1=aux_values[i_TOT - 2 - m][k][3][t]-aux_values[m][k][3][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=aux_values[m + 1][k][3][t]-aux_values[m][k][3][t];
							}
							if (k != kMax[0] - 2) {
								D_3=aux_values[m][k + 1][3][t]-aux_values[m][k][3][t];
							}
							if (m != 0) {
								D_4=aux_values[m - 1][k][3][t]-aux_values[m][k][3][t];
							}
							aux_values[m][k][7][t] = aux_values[m][k][7][t] + WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
							//---tau_zz
							D_1=0; D_2=0; D_3=0; D_4=0;	//Re-initialization
							if (k != 0) {
								D_1=aux_values[m][k - 1][6][t]-aux_values[m][k][6][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									//no action, D_1 equal zero (initialization) is ok
								} else {
									D_1=aux_values[i_TOT - 2 - m][k][6][t]-aux_values[m][k][6][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=aux_values[m + 1][k][6][t]-aux_values[m][k][6][t];
							}
							if (k != kMax[0] - 2) {
								D_3=aux_values[m][k + 1][6][t]-aux_values[m][k][6][t];
							}
							if (m != 0) {
								D_4=aux_values[m - 1][k][6][t]-aux_values[m][k][6][t];
							}
							aux_values[m][k][8][t] = aux_values[m][k][8][t] + WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
							break;
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of computing the divergence of the viscous stress tensor.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Transfer divergence of the viscous stress tensor from cells centres to grid nodes, or compute it at grid nodes in case of FD scheme selected.
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		for (int k = 0; k <kMax[0]; k++) {
			for (int m = 0; m <i_TOT; m++) {
				//TRANSFER from CC to GN
				if ((!flux_scheme)&&(grad_scheme!=1)) {		//Only transfer these variables in case they are needed (flux_scheme==0) and not directly computed at grid nodes (grad_scheme!=1)
					double inv_d1 = 0;
					double inv_d2 = 0;
					double inv_d3 = 0;
					double inv_d4 = 0;
					//---Divergence of the viscous stress tensor
					if (k==kMax[0]-1) {	//Node is on far-field boundary (C-edge)
						if (m==0) {
							aux_values_N[m][k][12][t] = aux_values[m][k - 1][7][t];
							aux_values_N[m][k][13][t] = aux_values[m][k - 1][8][t];

						} else if (m==i_TOT-1) {
							aux_values_N[m][k][12][t] = aux_values[m - 1][k - 1][7][t];
							aux_values_N[m][k][13][t] = aux_values[m - 1][k - 1][8][t];

						} else {
							inv_d1= 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d2= 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][12][t] = (aux_values[m - 1][k - 1][7][t] * inv_d1 + aux_values[m][k - 1][7][t] * inv_d2) / (inv_d1 + inv_d2);
							aux_values_N[m][k][13][t] = (aux_values[m - 1][k - 1][8][t] * inv_d1 + aux_values[m][k - 1][8][t] * inv_d2) / (inv_d1 + inv_d2);
						}
					} else if (k==0) {
						if (m==0) {
							inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[i_TOT - 2][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[i_TOT - 2][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][12][t] = (aux_values[i_TOT - 2][k][7][t] * inv_d2 + aux_values[m][k][7][t] * inv_d3) / (inv_d2 + inv_d3);
							aux_values_N[m][k][13][t] = (aux_values[i_TOT - 2][k][8][t] * inv_d2 + aux_values[m][k][8][t] * inv_d3) / (inv_d2 + inv_d3);
						} else if (m==i_TOT-1) {
							inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[0][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[0][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][12][t] = (aux_values[m - 1][k][7][t] * inv_d4 + aux_values[0][k][7][t] * inv_d1) / (inv_d1 + inv_d4);
							aux_values_N[m][k][13][t] = (aux_values[m - 1][k][8][t] * inv_d4 + aux_values[0][k][8][t] * inv_d1) / (inv_d1 + inv_d4);
						} else if (near_ID[m]) {	//Node is on body surface
							inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][12][t] = (aux_values[m - 1][k][7][t] * inv_d4 + aux_values[m][k][7][t] * inv_d3) / (inv_d3 + inv_d4);
							aux_values_N[m][k][13][t] = (aux_values[m - 1][k][8][t] * inv_d4 + aux_values[m][k][8][t] * inv_d3) / (inv_d3 + inv_d4);
						} else {	//Node is inside fluid domain
							inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[i_TOT - 1 - m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[i_TOT - 1 - m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
												pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
							aux_values_N[m][k][12][t] = (aux_values[i_TOT - 1 - m][k][7][t] * inv_d1 + aux_values[i_TOT - 1 - m - 1][k][7][t] * inv_d2 +
																aux_values[m][k][7][t] * inv_d3 + aux_values[m - 1][k][7][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
							aux_values_N[m][k][13][t] = (aux_values[i_TOT - 1 - m][k][8][t] * inv_d1 + aux_values[i_TOT - 1 - m - 1][k][8][t] * inv_d2 +
																aux_values[m][k][8][t] * inv_d3 + aux_values[m - 1][k][8][t] * inv_d4) /
																(inv_d1 + inv_d2 + inv_d3 + inv_d4);
						}
					} else if (m==0) {	//Node is on far-field boundary (outflow edge)
						inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						aux_values_N[m][k][12][t] = (aux_values[m][k - 1][7][t] * inv_d2 + aux_values[m][k][7][t] * inv_d3) / (inv_d2 + inv_d3);
						aux_values_N[m][k][13][t] = (aux_values[m][k - 1][8][t] * inv_d2 + aux_values[m][k][8][t] * inv_d3) / (inv_d2 + inv_d3);
					} else if (m==i_TOT-1) {
						inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						aux_values_N[m][k][12][t] = (aux_values[m - 1][k - 1][7][t] * inv_d1 + aux_values[m - 1][k][7][t] * inv_d4) / (inv_d1 + inv_d4);
						aux_values_N[m][k][13][t] = (aux_values[m - 1][k - 1][8][t] * inv_d1 + aux_values[m - 1][k][8][t] * inv_d4) / (inv_d1 + inv_d4);
					} else {	//Node is in the mesh interior (not on grid boundaries)
						inv_d1 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m - 1][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						inv_d2 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k - 1][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m][k - 1][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						inv_d3 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						inv_d4 = 1 / (sqrt(pow(nodal_values[m][k][0][t] - cell_values[m - 1][k][numVars + compVarsN + 15][t] - X_POLE[t], 2) +
											pow(nodal_values[m][k][2][t] - cell_values[m - 1][k][numVars + compVarsN + 16][t] - Z_POLE[t], 2)));
						aux_values_N[m][k][12][t] = (aux_values[m - 1][k - 1][7][t] * inv_d1 + aux_values[m][k - 1][7][t] * inv_d2 +
															aux_values[m][k][7][t] * inv_d3 + aux_values[m - 1][k][7][t] * inv_d4) /
															(inv_d1 + inv_d2 + inv_d3 + inv_d4);
						aux_values_N[m][k][13][t] = (aux_values[m - 1][k - 1][8][t] * inv_d1 + aux_values[m][k - 1][8][t] * inv_d2 +
															aux_values[m][k][8][t] * inv_d3 + aux_values[m - 1][k][8][t] * inv_d4) /
															(inv_d1 + inv_d2 + inv_d3 + inv_d4);
					}
				}
				//COMPUTE at CC
				if (grad_scheme==1) {
					//Local Jacobian
					double J11=0, J12=0, J21=0, J22=0, ddm=0, ddk=0, ddx_t11=0, ddx_t13=0, ddz_t13=0, ddz_t33=0;
					if (m==0) {
						J11=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m+1][k][0][t]-0.5*nodal_values[m+2][k][0][t];
						J12=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m+1][k][2][t]-0.5*nodal_values[m+2][k][2][t];
					} else if (m==i_TOT-1) {
						J11=0.5*nodal_values[m-2][k][0][t]-2.0*nodal_values[m-1][k][0][t]+1.5*nodal_values[m][k][0][t];
						J12=0.5*nodal_values[m-2][k][2][t]-2.0*nodal_values[m-1][k][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J11=0.5*nodal_values[m+1][k][0][t]-0.5*nodal_values[m-1][k][0][t];
                    	J12=0.5*nodal_values[m+1][k][2][t]-0.5*nodal_values[m-1][k][2][t];
					}
                    if (k==0) {
						J21=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k+2][0][t];
						J22=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k+2][2][t];
					} else if (k==kMax[0]-1) {
						J21=0.5*nodal_values[m][k-2][0][t]-2.0*nodal_values[m][k-1][0][t]+1.5*nodal_values[m][k][0][t];
						J22=0.5*nodal_values[m][k-2][2][t]-2.0*nodal_values[m][k-1][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J21=0.5*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k-1][0][t];
						J22=0.5*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k-1][2][t];
					}
                    double detJ=J11*J22-J12*J21;
					//tau_11 index gradients
					if (m==0) {
						ddm=-1.5*aux_values_N[m][k][9][t]+2.0*aux_values_N[m+1][k][9][t]-0.5*aux_values_N[m+2][k][9][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*aux_values_N[m-2][k][9][t]-2.0*aux_values_N[m-1][k][9][t]+1.5*aux_values_N[m][k][9][t];
					} else {
						ddm=0.5*aux_values_N[m+1][k][9][t]-0.5*aux_values_N[m-1][k][9][t];
					}
					if (k==0) {
						ddk=-1.5*aux_values_N[m][k][9][t]+2.0*aux_values_N[m][k+1][9][t]-0.5*aux_values_N[m][k+2][9][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*aux_values_N[m][k-2][9][t]-2.0*aux_values_N[m][k-1][9][t]+1.5*aux_values_N[m][k][9][t];
					} else {
						ddk=0.5*aux_values_N[m][k+1][9][t]-0.5*aux_values_N[m][k-1][9][t];
					}
					//tau_11 spatial gradients
					ddx_t11 = (J22*ddm - J12*ddk) / detJ;
					//tau_13 index gradients
					if (m==0) {
						ddm=-1.5*aux_values_N[m][k][10][t]+2.0*aux_values_N[m+1][k][10][t]-0.5*aux_values_N[m+2][k][10][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*aux_values_N[m-2][k][10][t]-2.0*aux_values_N[m-1][k][10][t]+1.5*aux_values_N[m][k][10][t];
					} else {
						ddm=0.5*aux_values_N[m+1][k][10][t]-0.5*aux_values_N[m-1][k][10][t];
					}
					if (k==0) {
						ddk=-1.5*aux_values_N[m][k][10][t]+2.0*aux_values_N[m][k+1][10][t]-0.5*aux_values_N[m][k+2][10][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*aux_values_N[m][k-2][10][t]-2.0*aux_values_N[m][k-1][10][t]+1.5*aux_values_N[m][k][10][t];
					} else {
						ddk=0.5*aux_values_N[m][k+1][10][t]-0.5*aux_values_N[m][k-1][10][t];
					}
					//tau_13 spatial gradients
					ddx_t13 = (J22*ddm - J12*ddk) / detJ;
					ddz_t13 = (-J21*ddm + J11*ddk) / detJ;
					//tau_33 index gradients
					if (m==0) {
						ddm=-1.5*aux_values_N[m][k][11][t]+2.0*aux_values_N[m+1][k][11][t]-0.5*aux_values_N[m+2][k][11][t];
					} else if (m==i_TOT-1) {
						ddm=0.5*aux_values_N[m-2][k][11][t]-2.0*aux_values_N[m-1][k][11][t]+1.5*aux_values_N[m][k][11][t];
					} else {
						ddm=0.5*aux_values_N[m+1][k][11][t]-0.5*aux_values_N[m-1][k][11][t];
					}
					if (k==0) {
						ddk=-1.5*aux_values_N[m][k][11][t]+2.0*aux_values_N[m][k+1][11][t]-0.5*aux_values_N[m][k+2][11][t];
					} else if (k==kMax[0]-1) {
						ddk=0.5*aux_values_N[m][k-2][11][t]-2.0*aux_values_N[m][k-1][11][t]+1.5*aux_values_N[m][k][11][t];
					} else {
						ddk=0.5*aux_values_N[m][k+1][11][t]-0.5*aux_values_N[m][k-1][11][t];
					}
					//tau_33 spatial gradients
					ddz_t33 = (-J21*ddm + J11*ddk) / detJ;
					//Viscous stress tensor divergence
					aux_values_N[m][k][12][t] = ddx_t11 + ddz_t13;
					aux_values_N[m][k][13][t] = ddx_t13 + ddz_t33;
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of transfering the divergence of the viscous stress tensor from cells centres to grid nodes.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Average the divergence of the viscous stress tensor at cells centres, in the case FD gradient method is selected.
	start = chrono::steady_clock::now();
	if ((grad_scheme==1)&&(Lamb_form>0)) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0] - 1; k++) {
				for (int m = 0; m < i_TOT - 1; m++) {
					//Specific total enthalpy per unit mass
					aux_values[m][k][7][t] =
						(aux_values_N[m][k][12][t] + aux_values_N[m + 1][k][12][t] + aux_values_N[m][k + 1][12][t] + aux_values_N[m + 1][k + 1][12][t]) / 4.0;
					//Specific entropy per unit mass
					aux_values[m][k][8][t] =
						(aux_values_N[m][k][13][t] + aux_values_N[m + 1][k][13][t] + aux_values_N[m][k + 1][13][t] + aux_values_N[m + 1][k + 1][13][t]) / 4.0;
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of averaging divergence of the viscous stress tensor at cells centres.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Re-definition of Lamb vector (and related auxiliary variables) based on Crocco's equation, in case LAMB_FORM > 0
	//---Vorticity is also re-defined and derived via inverse cross-product with fluid velocity
	start = chrono::steady_clock::now();
	if (Lamb_form>0) {
		if (grad_scheme==0) {
			for (int t=0; t<Nt; t++) {
				for (int k = 0; k <kMax[0]-1; k++) {
					for (int m = 0; m <i_TOT-1; m++) {
						cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
									rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
						//Entropy variation gradient
						double ds_dx = 0, ds_dz = 0;
						ds_dx=(R_air/cell_values[m][k][numVars+compVarsN][t])*
								(((aux_values_N[m+1][k][3][t]+aux_values_N[m][k][3][t])/2)*dz_1+
								((aux_values_N[m+1][k+1][3][t]+aux_values_N[m+1][k][3][t])/2)*dz_2+
								((aux_values_N[m][k+1][3][t]+aux_values_N[m+1][k+1][3][t])/2)*dz_3+
								((aux_values_N[m][k][3][t]+aux_values_N[m][k+1][3][t])/2)*dz_4);
						ds_dz=(R_air/cell_values[m][k][numVars+compVarsN][t])*
								((-(aux_values_N[m+1][k][3][t]+aux_values_N[m][k][3][t])/2)*dx_1+
								(-(aux_values_N[m+1][k+1][3][t]+aux_values_N[m+1][k][3][t])/2)*dx_2+
								(-(aux_values_N[m][k+1][3][t]+aux_values_N[m+1][k+1][3][t])/2)*dx_3+
								(-(aux_values_N[m][k][3][t]+aux_values_N[m][k+1][3][t])/2)*dx_4);
						//Total enthalpy gradient
						double dH_dx = 0, dH_dz = 0;
						dH_dx=(pow(V_inf, 2)/cell_values[m][k][numVars+compVarsN][t])*
								(((aux_values_N[m+1][k][4][t]+aux_values_N[m][k][4][t])/2)*dz_1+
								((aux_values_N[m+1][k+1][4][t]+aux_values_N[m+1][k][4][t])/2)*dz_2+
								((aux_values_N[m][k+1][4][t]+aux_values_N[m+1][k+1][4][t])/2)*dz_3+
								((aux_values_N[m][k][4][t]+aux_values_N[m][k+1][4][t])/2)*dz_4);
						dH_dz=(pow(V_inf, 2)/cell_values[m][k][numVars+compVarsN][t])*
								((-(aux_values_N[m+1][k][4][t]+aux_values_N[m][k][4][t])/2)*dx_1+
								(-(aux_values_N[m+1][k+1][4][t]+aux_values_N[m+1][k][4][t])/2)*dx_2+
								(-(aux_values_N[m][k+1][4][t]+aux_values_N[m+1][k+1][4][t])/2)*dx_3+
								(-(aux_values_N[m][k][4][t]+aux_values_N[m][k+1][4][t])/2)*dx_4);
						//Lamb vector
						if (U_FORM==0) {
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx + aux_values[m][k][7][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t]+aux_values_N[m+1][k+1][27][t]+aux_values_N[m][k+1][27][t]);
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz + aux_values[m][k][8][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t]+aux_values_N[m+1][k+1][28][t]+aux_values_N[m][k+1][28][t]);
						} else if (U_FORM==1) {
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx + aux_values[m][k][7][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][25][t]+aux_values_N[m+1][k][25][t]+aux_values_N[m+1][k+1][25][t]+aux_values_N[m][k+1][25][t]) +
								0.25*(aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t]+aux_values_N[m+1][k+1][27][t]+aux_values_N[m][k+1][27][t]);
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz + aux_values[m][k][8][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][26][t]+aux_values_N[m+1][k][26][t]+aux_values_N[m+1][k+1][26][t]+aux_values_N[m][k+1][26][t]) +
								0.25*(aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t]+aux_values_N[m+1][k+1][28][t]+aux_values_N[m][k+1][28][t]);
						} else if (U_FORM==2) {
							
						}
						//Transverse vorticity by cross-product inverse
						cell_values[m][k][numVars+compVarsN+12][t]=(-cell_values[m][k][vx_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t]+
																	 cell_values[m][k][vz_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t])/
																	 (2*cell_values[m][k][numVars+1][t]);
						if (Lamb_form==2) {	//In this case the Lamb vector (and, then, associated auxiliary variables) are computed using the negative-G vector field, instead of the full RHS of Crocco's equation
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx;
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz;
						}
						//Auxiliar variables (for easier elementary aerodynamic force calculation: FS_omega terms)
						aux_values[m][k][0][t]=(cell_values[m][k][numVars+compVarsN+15][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t]+
							cell_values[m][k][numVars+compVarsN+16][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t]);
						aux_values[m][k][1][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t];
						aux_values[m][k][2][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t];
					}
				}
			}
		} else if (grad_scheme==1) {	//Re-define on nodes
			//Nodel gradients of mass-specific total enthalpy and entropy
			for (int t=0; t<Nt; t++) {
				for (int k = 0; k <kMax[0]; k++) {
					for (int m = 0; m <i_TOT; m++) {
						//Local Jacobian
						double J11=0, J12=0, J21=0, J22=0, ddm=0, ddk=0;
						if (m==0) {
							J11=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m+1][k][0][t]-0.5*nodal_values[m+2][k][0][t];
							J12=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m+1][k][2][t]-0.5*nodal_values[m+2][k][2][t];
						} else if (m==i_TOT-1) {
							J11=0.5*nodal_values[m-2][k][0][t]-2.0*nodal_values[m-1][k][0][t]+1.5*nodal_values[m][k][0][t];
							J12=0.5*nodal_values[m-2][k][2][t]-2.0*nodal_values[m-1][k][2][t]+1.5*nodal_values[m][k][2][t];
						} else {
							J11=0.5*nodal_values[m+1][k][0][t]-0.5*nodal_values[m-1][k][0][t];
							J12=0.5*nodal_values[m+1][k][2][t]-0.5*nodal_values[m-1][k][2][t];
						}
						if (k==0) {
							J21=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k+2][0][t];
							J22=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k+2][2][t];
						} else if (k==kMax[0]-1) {
							J21=0.5*nodal_values[m][k-2][0][t]-2.0*nodal_values[m][k-1][0][t]+1.5*nodal_values[m][k][0][t];
							J22=0.5*nodal_values[m][k-2][2][t]-2.0*nodal_values[m][k-1][2][t]+1.5*nodal_values[m][k][2][t];
						} else {
							J21=0.5*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k-1][0][t];
							J22=0.5*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k-1][2][t];
						}
						double detJ=J11*J22-J12*J21;
						//Specific total enthalpy index gradients
						if (m==0) {
							ddm=(-1.5*aux_values_N[m][k][4][t]+2.0*aux_values_N[m+1][k][4][t]-0.5*aux_values_N[m+2][k][4][t])*pow(V_inf, 2);
						} else if (m==i_TOT-1) {
							ddm=(0.5*aux_values_N[m-2][k][4][t]-2.0*aux_values_N[m-1][k][4][t]+1.5*aux_values_N[m][k][4][t])*pow(V_inf, 2);
						} else {
							ddm=(0.5*aux_values_N[m+1][k][4][t]-0.5*aux_values_N[m-1][k][4][t])*pow(V_inf, 2);
						}
						if (k==0) {
							ddk=(-1.5*aux_values_N[m][k][4][t]+2.0*aux_values_N[m][k+1][4][t]-0.5*aux_values_N[m][k+2][4][t])*pow(V_inf, 2);
						} else if (k==kMax[0]-1) {
							ddk=(0.5*aux_values_N[m][k-2][4][t]-2.0*aux_values_N[m][k-1][4][t]+1.5*aux_values_N[m][k][4][t])*pow(V_inf, 2);
						} else {
							ddk=(0.5*aux_values_N[m][k+1][4][t]-0.5*aux_values_N[m][k-1][4][t])*pow(V_inf, 2);
						}
						//Specific total enthalpy spatial gradients
						aux_values_N[m][k][32][t] = (J22*ddm - J12*ddk) / detJ;
						aux_values_N[m][k][33][t] = (-J21*ddm + J11*ddk) / detJ;
						//Specific entropy index gradients
						if (m==0) {
							ddm=(-1.5*aux_values_N[m][k][3][t]+2.0*aux_values_N[m+1][k][3][t]-0.5*aux_values_N[m+2][k][3][t])*R_air;
						} else if (m==i_TOT-1) {
							ddm=(0.5*aux_values_N[m-2][k][3][t]-2.0*aux_values_N[m-1][k][3][t]+1.5*aux_values_N[m][k][3][t])*R_air;
						} else {
							ddm=(0.5*aux_values_N[m+1][k][3][t]-0.5*aux_values_N[m-1][k][3][t])*R_air;
						}
						if (k==0) {
							ddk=(-1.5*aux_values_N[m][k][3][t]+2.0*aux_values_N[m][k+1][3][t]-0.5*aux_values_N[m][k+2][3][t])*R_air;
						} else if (k==kMax[0]-1) {
							ddk=(0.5*aux_values_N[m][k-2][3][t]-2.0*aux_values_N[m][k-1][3][t]+1.5*aux_values_N[m][k][3][t])*R_air;
						} else {
							ddk=(0.5*aux_values_N[m][k+1][3][t]-0.5*aux_values_N[m][k-1][3][t])*R_air;
						}
						//Specific entropy spatial gradients
						aux_values_N[m][k][34][t] = (J22*ddm - J12*ddk) / detJ;
						aux_values_N[m][k][35][t] = (-J21*ddm + J11*ddk) / detJ;
					}
				}
			}
			//---Calculate Lamb vector (and related auxiliary variables) based on Crocco's equation, using average values of entropy and total enthalpy gradients from grid nodes to cell-centres.
			for (int t=0; t<Nt; t++) {
				for (int k = 0; k < kMax[0] - 1; k++) {
					for (int m = 0; m < i_TOT - 1; m++) {
						double dH_dx = 0, dH_dz = 0, ds_dx = 0, ds_dz = 0;
						//Averaged gradient of specific total enthalpy per unit mass
						dH_dx =	(aux_values_N[m][k][32][t] + aux_values_N[m + 1][k][32][t] + aux_values_N[m][k + 1][32][t] + aux_values_N[m + 1][k + 1][32][t]) / 4.0;
						dH_dz =	(aux_values_N[m][k][33][t] + aux_values_N[m + 1][k][33][t] + aux_values_N[m][k + 1][33][t] + aux_values_N[m + 1][k + 1][33][t]) / 4.0;
						//averaged gradient of specific entropy per unit mass
						ds_dx =	(aux_values_N[m][k][34][t] + aux_values_N[m + 1][k][34][t] + aux_values_N[m][k + 1][34][t] + aux_values_N[m + 1][k + 1][34][t]) / 4.0;
						ds_dz =	(aux_values_N[m][k][35][t] + aux_values_N[m + 1][k][35][t] + aux_values_N[m][k + 1][35][t] + aux_values_N[m + 1][k + 1][35][t]) / 4.0;
						//Lamb vector
						if (U_FORM==0) {
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx + aux_values[m][k][7][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t]+aux_values_N[m+1][k+1][27][t]+aux_values_N[m][k+1][27][t]);
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz + aux_values[m][k][8][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t]+aux_values_N[m+1][k+1][28][t]+aux_values_N[m][k+1][28][t]);
						} else if (U_FORM==1) {
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx + aux_values[m][k][7][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][25][t]+aux_values_N[m+1][k][25][t]+aux_values_N[m+1][k+1][25][t]+aux_values_N[m][k+1][25][t]) +
								0.25*(aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t]+aux_values_N[m+1][k+1][27][t]+aux_values_N[m][k+1][27][t]);
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz + aux_values[m][k][8][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][26][t]+aux_values_N[m+1][k][26][t]+aux_values_N[m+1][k+1][26][t]+aux_values_N[m][k+1][26][t]) +
								0.25*(aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t]+aux_values_N[m+1][k+1][28][t]+aux_values_N[m][k+1][28][t]);
						} else if (U_FORM==2) {
							
						}
						//Transverse vorticity by cross-product inverse
						cell_values[m][k][numVars+compVarsN+12][t]=(-cell_values[m][k][vx_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t]+
																	 cell_values[m][k][vz_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t])/
																	 (2*cell_values[m][k][numVars+1][t]);
						if (Lamb_form==2) {	//In this case the Lamb vector (and, then, associated auxiliary variables) are computed using the negative-G vector field, instead of the full RHS of Crocco's equation
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx;
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz;
						}
						//Auxiliar variables (for easier elementary aerodynamic force calculation: FS_omega terms)
						aux_values[m][k][0][t]=(cell_values[m][k][numVars+compVarsN+15][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t]+
							cell_values[m][k][numVars+compVarsN+16][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t]);
						aux_values[m][k][1][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t];
						aux_values[m][k][2][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t];
					}
				}
			}
		} else if (grad_scheme==2) {
			for (int t=0; t<Nt; t++) {
				for (int k = 0; k <kMax[0]-1; k++) {
					for (int m = 0; m <i_TOT-1; m++) {
						cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
									rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
						double D_1, D_2, D_3, D_4;
						double dH_dx = 0, dH_dz = 0, ds_dx = 0, ds_dz = 0;
						//---Specific total enthalpy per unit mass
							D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
							if (k != 0) {
								D_1=aux_values[m][k - 1][14][t]-aux_values[m][k][14][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									D_1=-(aux_values[m][k + 1][14][t]-aux_values[m][k][14][t]);
								} else {
									D_1=aux_values[i_TOT - 2 - m][k][14][t]-aux_values[m][k][14][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=aux_values[m + 1][k][14][t]-aux_values[m][k][14][t];
							}
							if (k != kMax[0] - 2) {
								D_3=aux_values[m][k + 1][14][t]-aux_values[m][k][14][t];
							}
							if (m != 0) {
								D_4=aux_values[m - 1][k][14][t]-aux_values[m][k][14][t];
							}
						dH_dx=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
						dH_dz=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						//---Specific entropy per unit mass
							D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
							if (k != 0) {
								D_1=aux_values[m][k - 1][15][t]-aux_values[m][k][15][t];
							} else {
								if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
									D_1=-(aux_values[m][k + 1][15][t]-aux_values[m][k][15][t]);
								} else {
									D_1=aux_values[i_TOT - 2 - m][k][15][t]-aux_values[m][k][15][t];
								}
							}
							if (m != i_TOT - 2) {
								D_2=aux_values[m + 1][k][15][t]-aux_values[m][k][15][t];
							}
							if (k != kMax[0] - 2) {
								D_3=aux_values[m][k + 1][15][t]-aux_values[m][k][15][t];
							}
							if (m != 0) {
								D_4=aux_values[m - 1][k][15][t]-aux_values[m][k][15][t];
							}
						ds_dx=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
						ds_dz=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
						//Lamb vector
						if (U_FORM==0) {
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx + aux_values[m][k][7][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t]+aux_values_N[m+1][k+1][27][t]+aux_values_N[m][k+1][27][t]);
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz + aux_values[m][k][8][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t]+aux_values_N[m+1][k+1][28][t]+aux_values_N[m][k+1][28][t]);
						} else if (U_FORM==1) {
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx + aux_values[m][k][7][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][25][t]+aux_values_N[m+1][k][25][t]+aux_values_N[m+1][k+1][25][t]+aux_values_N[m][k+1][25][t]) +
								0.25*(aux_values_N[m][k][27][t]+aux_values_N[m+1][k][27][t]+aux_values_N[m+1][k+1][27][t]+aux_values_N[m][k+1][27][t]);
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz + aux_values[m][k][8][t]/cell_values[m][k][rho_id-1][t] -
								0.25*(aux_values_N[m][k][26][t]+aux_values_N[m+1][k][26][t]+aux_values_N[m+1][k+1][26][t]+aux_values_N[m][k+1][26][t]) +
								0.25*(aux_values_N[m][k][28][t]+aux_values_N[m+1][k][28][t]+aux_values_N[m+1][k+1][28][t]+aux_values_N[m][k+1][28][t]);
						} else if (U_FORM==2) {
							
						}
						//Transverse vorticity by cross-product inverse
						cell_values[m][k][numVars+compVarsN+12][t]=(-cell_values[m][k][vx_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t]+
																	 cell_values[m][k][vz_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t])/
																	 (2*cell_values[m][k][numVars+1][t]);
						if (Lamb_form==2) {	//In this case the Lamb vector (and, then, associated auxiliary variables) are computed using the negative-G vector field, instead of the full RHS of Crocco's equation
							cell_values[m][k][numVars+compVarsN+13][t] = cell_values[m][k][T_id-1][t]*ds_dx - dH_dx;
							cell_values[m][k][numVars+compVarsN+14][t] = cell_values[m][k][T_id-1][t]*ds_dz - dH_dz;
						}
						//Auxiliar variables (for easier elementary aerodynamic force calculation: FS_omega terms)
						aux_values[m][k][0][t]=(cell_values[m][k][numVars+compVarsN+15][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t]+
							cell_values[m][k][numVars+compVarsN+16][t]*cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t]);
						aux_values[m][k][1][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t];
						aux_values[m][k][2][t]=cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t];
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of Lamb vector re-definition based on Crocco's equation.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Average values of the divergence of the viscous stress tensor at cells centres, in the case FD gradient method is selected.
	start = chrono::steady_clock::now();
	if (grad_scheme==1) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0] - 1; k++) {
				for (int m = 0; m < i_TOT - 1; m++) {
				aux_values[m][k][7][t] =
					(aux_values_N[m][k][12][t] + aux_values_N[m + 1][k][12][t] + aux_values_N[m][k + 1][12][t] + aux_values_N[m + 1][k + 1][12][t]) / 4.0;
				aux_values[m][k][8][t] =
					(aux_values_N[m][k][13][t] + aux_values_N[m + 1][k][13][t] + aux_values_N[m][k + 1][13][t] + aux_values_N[m + 1][k + 1][13][t]) / 4.0;
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of averaging values the divergence of the viscous stress tensor at cells centres.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Average values of time-derivatives nodal variables at cells centres, in the case WLS gradient method is selected and FT_FORM==false.
	//---TEMPORARY auxiliary variables are used (and reused for other quantities afterwards in the code...)
	start = chrono::steady_clock::now();
	if ((grad_scheme==2)&&(!FT_FORM)) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0] - 1; k++) {
				for (int m = 0; m < i_TOT - 1; m++) {
					if (U_FORM==0) {				//Average values of (rho*dV/dt)
						aux_values[m][k][14][t] =
							(aux_values_N[m][k][27][t] + aux_values_N[m + 1][k][27][t] + aux_values_N[m][k + 1][27][t] + aux_values_N[m + 1][k + 1][27][t]) / 4.0;
						aux_values[m][k][15][t] =
							(aux_values_N[m][k][28][t] + aux_values_N[m + 1][k][28][t] + aux_values_N[m][k + 1][28][t] + aux_values_N[m + 1][k + 1][28][t]) / 4.0;
					} else if ((U_FORM==1)||(U_FORM==2)) {			//Average values of (rho*dV'/dt') or (rho*dV/dt')
						aux_values[m][k][14][t] =
							(aux_values_N[m][k][25][t] + aux_values_N[m + 1][k][25][t] + aux_values_N[m][k + 1][25][t] + aux_values_N[m + 1][k + 1][25][t]) / 4.0;
						aux_values[m][k][15][t] =
							(aux_values_N[m][k][26][t] + aux_values_N[m + 1][k][26][t] + aux_values_N[m][k + 1][26][t] + aux_values_N[m + 1][k + 1][26][t]) / 4.0;
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of averaging values of time-derivatives nodal variables.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}
	//---AERODYNAMIC FORCE CALCULATION
	start = chrono::steady_clock::now();
	vector < double > Clp_nf(Nt,0);
	vector < double > Clf_nf(Nt,0);
	vector < double > Cdp_nf(Nt,0);
	vector < double > Cdf_nf(Nt,0);
	for (int t=0; t<Nt; t++) {
		double CdeP_thd_1 = 0;
		double CdeP_thd_2 = 0;
		for (int k = 0; k <kMax[0]-1; k++) {
			for (int m = 0; m <i_TOT-1; m++) {
				//---Metrics of elementary cell surfaces
				cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
							rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
				//Check values of fluxes blend coefficients
				if ((a_1 < 0) || (a_1 > 1) || (a_2 < 0) || (a_2 > 1) || (a_3 < 0) || (a_3 > 1) || (a_4 < 0) || (a_4 > 1)) {
					cout << m << " " << k << " " << a_1 << " " << a_2 << " " << a_3 << " " << a_4 << "\n";
				}
				//NEAR-FIELD aerodynamic force re-computation
				if ((k == 0) && ((near_ID[m]) || (near_ID[m + 1]))) {
					double T_fc = (nodal_values[m][k][T_id - 1][t] + nodal_values[m + 1][k][T_id - 1][t]) / 2;
					Clp_nf[t] = Clp_nf[t] - ((nodal_values[m][k][numVars][t] + nodal_values[m + 1][k][numVars][t]) / 2) * dx_1;
					Cdp_nf[t] = Cdp_nf[t] + ((nodal_values[m][k][numVars][t] + nodal_values[m + 1][k][numVars][t]) / 2) * dz_1;
					if (grad_scheme!=1) {		//Extrapolate gradient from CC of first cell adjacent to the wall (a test with two-cells gradient extrapolation was done and commented-out)
						Clf_nf[t] = Clf_nf[t] - 2 * (((nodal_values[m][k][mut_id - 1][t] + nodal_values[m + 1][k][mut_id - 1][t]) / 2) + ((0.000001458 * sqrt(T_fc)) / (1 + 110.4 / T_fc))) *
							(((cell_values[m][k][numVars + compVarsN + 3][t] + cell_values[m][k][numVars + compVarsN + 4][t]) / 2) * dz_1 - dx_1 *
								(cell_values[m][k][numVars + compVarsN + 5][t] - (1.0 / 3.0) * (cell_values[m][k][numVars + compVarsN + 2][t] + cell_values[m][k][numVars + compVarsN + 5][t])));
						Cdf_nf[t] = Cdf_nf[t] - 2 * (((nodal_values[m][k][mut_id - 1][t] + nodal_values[m + 1][k][mut_id - 1][t]) / 2) + ((0.000001458 * sqrt(T_fc)) / (1 + 110.4 / T_fc))) *
							((cell_values[m][k][numVars + compVarsN + 2][t] - (1.0 / 3.0) * (cell_values[m][k][numVars + compVarsN + 2][t] + cell_values[m][k][numVars + compVarsN + 5][t])) * dz_1 -
								dx_1 * ((cell_values[m][k][numVars + compVarsN + 3][t] + cell_values[m][k][numVars + compVarsN + 4][t]) / 2));
						//Cdf_nf = Cdf_nf - 
						//	((-a_1 * aux_values[m][k+1][3] + (1 + a_1) * aux_values[m][k][3]) * dz_1 -
						//		dx_1 * (-a_1 * aux_values[m][k+1][4] + (1 + a_1) * aux_values[m][k][4]));
					} else {	//Use velocity gradients at wall grid nodes
						Clf_nf[t] = Clf_nf[t] - 2 * (((nodal_values[m][k][mut_id - 1][t] + nodal_values[m + 1][k][mut_id - 1][t]) / 2) + ((0.000001458 * sqrt(T_fc)) / (1 + 110.4 / T_fc))) *
							((((aux_values_N[m][k][22][t]+aux_values_N[m+1][k][22][t])*0.5 + (aux_values_N[m][k][23][t]+aux_values_N[m+1][k][23][t])*0.5) / 2) * dz_1 -
							dx_1 * ((aux_values_N[m][k][24][t]+aux_values_N[m+1][k][24][t])*0.5 - (1.0 / 3.0) *
								   ((aux_values_N[m][k][21][t]+aux_values_N[m+1][k][21][t])*0.5 + (aux_values_N[m][k][24][t]+aux_values_N[m+1][k][24][t])*0.5)));
						Cdf_nf[t] = Cdf_nf[t] - 2 * (((nodal_values[m][k][mut_id - 1][t] + nodal_values[m + 1][k][mut_id - 1][t]) / 2) + ((0.000001458 * sqrt(T_fc)) / (1 + 110.4 / T_fc))) *
							(((aux_values_N[m][k][21][t]+aux_values_N[m+1][k][21][t])*0.5 - (1.0 / 3.0) *
							((aux_values_N[m][k][21][t]+aux_values_N[m+1][k][21][t])*0.5 + (aux_values_N[m][k][24][t]+aux_values_N[m+1][k][24][t])*0.5)) * dz_1 -
								dx_1 * (((aux_values_N[m][k][22][t]+aux_values_N[m+1][k][22][t])*0.5 + (aux_values_N[m][k][23][t]+aux_values_N[m+1][k][23][t])*0.5) / 2));
					}
				}
				//FAR-FIELD elementary aerodynamic force using vorticity-based methods
				//---Unsteady term (volume integrand)
				if (!FT_FORM) {
					//------Spatial gradients of (rho*dV/dt) or (rho*dV'/dt') or (rho*dV/dt')
					double drhoVtx_dz=NAN, drhoVtz_dx=NAN;
					if (U_FORM==0) {				//Spatial gradients of (rho*dV/dt)
						switch (grad_scheme) {
							case 0:
								drhoVtx_dz=(1/cell_values[m][k][numVars+compVarsN][t])*
										((-(nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][27][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][27][t])/2)*dx_1+
										(-(nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][27][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][27][t])/2)*dx_2+
										(-(nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][27][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][27][t])/2)*dx_3+
										(-(nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][27][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][27][t])/2)*dx_4);
								drhoVtz_dx=(1/cell_values[m][k][numVars+compVarsN][t])*
										(((nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][28][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][28][t])/2)*dz_1+
										((nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][28][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][28][t])/2)*dz_2+
										((nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][28][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][28][t])/2)*dz_3+
										((nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][28][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][28][t])/2)*dz_4);
								break;
							case 1:		//Do nothing. spatial gradients will be calculated at grid nodes later on, then averaged to CC in a following step.
								break;
							case 2:		//Weighted-Least-Squares
								double D_1, D_2, D_3, D_4;
								//---drhoVtx_dz
								D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=aux_values[m][k - 1][14][t]-aux_values[m][k][14][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*aux_values[m][k][14][t];
									} else {
										D_1=aux_values[i_TOT - 2 - m][k][14][t]-aux_values[m][k][14][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=aux_values[m + 1][k][14][t]-aux_values[m][k][14][t];
								}
								if (k != kMax[0] - 2) {
									D_3=aux_values[m][k + 1][14][t]-aux_values[m][k][14][t];
								}
								if (m != 0) {
									D_4=aux_values[m - 1][k][14][t]-aux_values[m][k][14][t];
								}
								drhoVtx_dz=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
								//---drhoVtz_dx
								D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=aux_values[m][k - 1][15][t]-aux_values[m][k][15][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*aux_values[m][k][15][t];
									} else {
										D_1=aux_values[i_TOT - 2 - m][k][15][t]-aux_values[m][k][15][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=aux_values[m + 1][k][15][t]-aux_values[m][k][15][t];
								}
								if (k != kMax[0] - 2) {
									D_3=aux_values[m][k + 1][15][t]-aux_values[m][k][15][t];
								}
								if (m != 0) {
									D_4=aux_values[m - 1][k][15][t]-aux_values[m][k][15][t];
								}
								drhoVtz_dx=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
								break;
						}
					} else if (U_FORM==1) {			//Spatial gradients of (rho*dV'/dt')
						switch (grad_scheme) {
							case 0:
								drhoVtx_dz=(1/cell_values[m][k][numVars+compVarsN][t])*
										((-(nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][25][t])/2)*dx_1+
										(-(nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][25][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t])/2)*dx_2+
										(-(nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][25][t])/2)*dx_3+
										(-(nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][25][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t])/2)*dx_4);
								drhoVtz_dx=(1/cell_values[m][k][numVars+compVarsN][t])*
										(((nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][26][t])/2)*dz_1+
										((nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][26][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t])/2)*dz_2+
										((nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][26][t])/2)*dz_3+
										((nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][26][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t])/2)*dz_4);
								break;
							case 1:		//Do nothing. spatial gradients will be calculated at grid nodes later on, then averaged to CC in a following step.
								break;
							case 2:		//Weighted-Least-Squares
								double D_1, D_2, D_3, D_4;
								//---drhoVtx_dz
								D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=aux_values[m][k - 1][14][t]-aux_values[m][k][14][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*aux_values[m][k][14][t];
									} else {
										D_1=aux_values[i_TOT - 2 - m][k][14][t]-aux_values[m][k][14][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=aux_values[m + 1][k][14][t]-aux_values[m][k][14][t];
								}
								if (k != kMax[0] - 2) {
									D_3=aux_values[m][k + 1][14][t]-aux_values[m][k][14][t];
								}
								if (m != 0) {
									D_4=aux_values[m - 1][k][14][t]-aux_values[m][k][14][t];
								}
								drhoVtx_dz=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
								//---drhoVtz_dx
								D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=aux_values[m][k - 1][15][t]-aux_values[m][k][15][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*aux_values[m][k][15][t];
									} else {
										D_1=aux_values[i_TOT - 2 - m][k][15][t]-aux_values[m][k][15][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=aux_values[m + 1][k][15][t]-aux_values[m][k][15][t];
								}
								if (k != kMax[0] - 2) {
									D_3=aux_values[m][k + 1][15][t]-aux_values[m][k][15][t];
								}
								if (m != 0) {
									D_4=aux_values[m - 1][k][15][t]-aux_values[m][k][15][t];
								}
								drhoVtz_dx=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
								break;
						}
					} else if (U_FORM==2) {			//Spatial gradients of (rho*dV/dt')
						switch (grad_scheme) {
							case 0:
								drhoVtx_dz=(1/cell_values[m][k][numVars+compVarsN][t])*
										((-(nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][25][t])/2)*dx_1+
										(-(nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][25][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t])/2)*dx_2+
										(-(nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][25][t])/2)*dx_3+
										(-(nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][25][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t])/2)*dx_4);
								drhoVtz_dx=(1/cell_values[m][k][numVars+compVarsN][t])*
										(((nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][26][t])/2)*dz_1+
										((nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][26][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t])/2)*dz_2+
										((nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][26][t])/2)*dz_3+
										((nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][26][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t])/2)*dz_4);
								break;
							case 1:		//Do nothing. spatial gradients will be calculated at grid nodes later on, then averaged to CC in a following step.
								break;
							case 2:		//Weighted-Least-Squares
								double D_1, D_2, D_3, D_4;
								//---drhoVtx_dz
								D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=aux_values[m][k - 1][14][t]-aux_values[m][k][14][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*aux_values[m][k][14][t];
									} else {
										D_1=aux_values[i_TOT - 2 - m][k][14][t]-aux_values[m][k][14][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=aux_values[m + 1][k][14][t]-aux_values[m][k][14][t];
								}
								if (k != kMax[0] - 2) {
									D_3=aux_values[m][k + 1][14][t]-aux_values[m][k][14][t];
								}
								if (m != 0) {
									D_4=aux_values[m - 1][k][14][t]-aux_values[m][k][14][t];
								}
								drhoVtx_dz=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
								//---drhoVtz_dx
								D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=aux_values[m][k - 1][15][t]-aux_values[m][k][15][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*aux_values[m][k][15][t];
									} else {
										D_1=aux_values[i_TOT - 2 - m][k][15][t]-aux_values[m][k][15][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=aux_values[m + 1][k][15][t]-aux_values[m][k][15][t];
								}
								if (k != kMax[0] - 2) {
									D_3=aux_values[m][k + 1][15][t]-aux_values[m][k][15][t];
								}
								if (m != 0) {
									D_4=aux_values[m - 1][k][15][t]-aux_values[m][k][15][t];
								}
								drhoVtz_dx=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
								break;
						}
					}
					//------Unsteady volume integrand
					if (grad_scheme!=1) {	//Save CPU as these variables will be effectively computed in a following code section
						AF_elem[m][k][35][t] = cell_values[m][k][numVars + compVarsN + 16][t] * (drhoVtx_dz - drhoVtz_dx) * cell_values[m][k][numVars + compVarsN][t];
						AF_elem[m][k][36][t] = -cell_values[m][k][numVars + compVarsN + 15][t] * (drhoVtx_dz - drhoVtz_dx) * cell_values[m][k][numVars + compVarsN][t];
					}
				} else {
					if (U_FORM==0) {		//Use inertial time derivative of absolute fluid velocity
						AF_elem[m][k][35][t] = - 0.25*(aux_values_N[m][k][27][t]*nodal_values[m][k][rho_id-1][t]+aux_values_N[m+1][k][27][t]*nodal_values[m+1][k][rho_id-1][t]+
											aux_values_N[m+1][k+1][27][t]*nodal_values[m+1][k+1][rho_id-1][t]+aux_values_N[m][k+1][27][t]*nodal_values[m][k+1][rho_id-1][t])*
											cell_values[m][k][numVars + compVarsN][t];
						AF_elem[m][k][36][t] = - 0.25*(aux_values_N[m][k][28][t]*nodal_values[m][k][rho_id-1][t]+aux_values_N[m+1][k][28][t]*nodal_values[m+1][k][rho_id-1][t]+
											aux_values_N[m+1][k+1][28][t]*nodal_values[m+1][k+1][rho_id-1][t]+aux_values_N[m][k+1][28][t]*nodal_values[m][k+1][rho_id-1][t])*
											cell_values[m][k][numVars + compVarsN][t];
					} else if (U_FORM==1) {		//Use non-inertial time derivative of relative fluid velocity
						AF_elem[m][k][35][t] = - 0.25*(aux_values_N[m][k][25][t]*nodal_values[m][k][rho_id-1][t]+aux_values_N[m+1][k][25][t]*nodal_values[m+1][k][rho_id-1][t]+
											aux_values_N[m+1][k+1][25][t]*nodal_values[m+1][k+1][rho_id-1][t]+aux_values_N[m][k+1][25][t]*nodal_values[m][k+1][rho_id-1][t])*
											cell_values[m][k][numVars + compVarsN][t];
						AF_elem[m][k][36][t] = - 0.25*(aux_values_N[m][k][26][t]*nodal_values[m][k][rho_id-1][t]+aux_values_N[m+1][k][26][t]*nodal_values[m+1][k][rho_id-1][t]+
											aux_values_N[m+1][k+1][26][t]*nodal_values[m+1][k+1][rho_id-1][t]+aux_values_N[m][k+1][26][t]*nodal_values[m][k+1][rho_id-1][t])*
											cell_values[m][k][numVars + compVarsN][t];
					} else if (U_FORM==2) {
						//To be implemented...
					}
				}
				//------Conditionally add turbulent term in the case of LAMB_FORM=1 or LAMB_FORM=2
				if (Lamb_form>0) {
				//------Spatial gradients of (rho*grad(kt))
					double drhoGradXkt_dz=NAN, drhoGradZkt_dx=NAN;
					switch (grad_scheme) {
						case 0:
							drhoGradXkt_dz=(1/cell_values[m][k][numVars+compVarsN][t])*
									((-(nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][30][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][30][t])/2)*dx_1+
									(-(nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][30][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][30][t])/2)*dx_2+
									(-(nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][30][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][30][t])/2)*dx_3+
									(-(nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][30][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][30][t])/2)*dx_4);
							drhoGradZkt_dx=(1/cell_values[m][k][numVars+compVarsN][t])*
									(((nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][31][t]+nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][31][t])/2)*dz_1+
									((nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][31][t]+nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][31][t])/2)*dz_2+
									((nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][31][t]+nodal_values[m+1][k+1][rho_id-1][t]*aux_values_N[m+1][k+1][31][t])/2)*dz_3+
									((nodal_values[m][k][rho_id-1][t]*aux_values_N[m][k][31][t]+nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][31][t])/2)*dz_4);
							break;
						case 1:		//Do nothing. spatial gradients will be calculated at grid nodes later on, then averaged to CC in a following step.
							break;
						case 2:		//Weighted-Least-Squares
							double D_1, D_2, D_3, D_4;
							//---drhoGradXkt_dz
							D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=cell_values[m][k - 1][rho_id-1][t]*aux_values[m][k - 1][12][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][12][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*cell_values[m][k][rho_id-1][t]*aux_values[m][k][12][t];
									} else {
										D_1=cell_values[i_TOT - 2 - m][k][rho_id-1][t]*aux_values[i_TOT - 2 - m][k][12][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][12][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=cell_values[m + 1][k][rho_id-1][t]*aux_values[m + 1][k][12][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][12][t];
								}
								if (k != kMax[0] - 2) {
									D_3=cell_values[m][k + 1][rho_id-1][t]*aux_values[m][k + 1][12][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][12][t];
								}
								if (m != 0) {
									D_4=cell_values[m - 1][k][rho_id-1][t]*aux_values[m - 1][k][12][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][12][t];
								}
								drhoGradXkt_dz=WLS[1][0]*D_1 + WLS[1][1]*D_2 + WLS[1][2]*D_3 + WLS[1][3]*D_4;
							//---drhoGradZkt_dx
							D_1=0; D_2=0; D_3=0; D_4=0;	//Initialization
								if (k != 0) {
									D_1=cell_values[m][k - 1][rho_id-1][t]*aux_values[m][k - 1][13][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][13][t];
								} else {
									if ((near_ID[m]) || (near_ID[m + 1])) {	//cell is adjacent to the body surface
										D_1=-2*cell_values[m][k][rho_id-1][t]*aux_values[m][k][13][t];
									} else {
										D_1=cell_values[i_TOT - 2 - m][k][rho_id-1][t]*aux_values[i_TOT - 2 - m][k][13][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][13][t];
									}
								}
								if (m != i_TOT - 2) {
									D_2=cell_values[m + 1][k][rho_id-1][t]*aux_values[m + 1][k][13][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][13][t];
								}
								if (k != kMax[0] - 2) {
									D_3=cell_values[m][k + 1][rho_id-1][t]*aux_values[m][k + 1][13][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][13][t];
								}
								if (m != 0) {
									D_4=cell_values[m - 1][k][rho_id-1][t]*aux_values[m - 1][k][13][t]-cell_values[m][k][rho_id-1][t]*aux_values[m][k][13][t];
								}
								drhoGradZkt_dx=WLS[0][0]*D_1 + WLS[0][1]*D_2 + WLS[0][2]*D_3 + WLS[0][3]*D_4;
								break;
						}
					//------Unsteady volume integrand
					if (grad_scheme!=1) {	//Save CPU as these variables will be effectively computed in a following code section
						AF_elem[m][k][35][t] = AF_elem[m][k][35][t] + cell_values[m][k][numVars + compVarsN + 16][t] * (drhoGradXkt_dz - drhoGradZkt_dx) * cell_values[m][k][numVars + compVarsN][t];
						AF_elem[m][k][36][t] = AF_elem[m][k][36][t] -cell_values[m][k][numVars + compVarsN + 15][t] * (drhoGradXkt_dz - drhoGradZkt_dx) * cell_values[m][k][numVars + compVarsN][t];
					}
				}
				//---Unsteady volume integrand of Formulation C
				//if (!INC) {
					AF_elem[m][k][37][t] = - 0.25*(aux_values_N[m][k][29][t]*nodal_values[m][k][vx_id-1][t]+aux_values_N[m+1][k][29][t]*nodal_values[m+1][k][vx_id-1][t]+
											aux_values_N[m+1][k+1][29][t]*nodal_values[m+1][k+1][vx_id-1][t]+aux_values_N[m][k+1][29][t]*nodal_values[m][k+1][vx_id-1][t])*
											cell_values[m][k][numVars + compVarsN][t];
					AF_elem[m][k][38][t] = - 0.25*(aux_values_N[m][k][29][t]*nodal_values[m][k][vz_id-1][t]+aux_values_N[m+1][k][29][t]*nodal_values[m+1][k][vz_id-1][t]+
											aux_values_N[m+1][k+1][29][t]*nodal_values[m+1][k+1][vz_id-1][t]+aux_values_N[m][k+1][29][t]*nodal_values[m][k+1][vz_id-1][t])*
											cell_values[m][k][numVars + compVarsN][t];
				//}
				//---Unsteady volume integrand of apparent force from relative motion term
				if (U_FORM==1) {
					AF_elem[m][k][39][t] = 0.25*(aux_values_N[m][k][27][t]*nodal_values[m][k][rho_id-1][t]+aux_values_N[m+1][k][27][t]*nodal_values[m+1][k][rho_id-1][t]+
												aux_values_N[m+1][k+1][27][t]*nodal_values[m+1][k+1][rho_id-1][t]+aux_values_N[m][k+1][27][t]*nodal_values[m][k+1][rho_id-1][t])*
												cell_values[m][k][numVars + compVarsN][t];
					AF_elem[m][k][40][t] = 0.25*(aux_values_N[m][k][28][t]*nodal_values[m][k][rho_id-1][t]+aux_values_N[m+1][k][28][t]*nodal_values[m+1][k][rho_id-1][t]+
												aux_values_N[m+1][k+1][28][t]*nodal_values[m+1][k+1][rho_id-1][t]+aux_values_N[m][k+1][28][t]*nodal_values[m][k+1][rho_id-1][t])*
												cell_values[m][k][numVars + compVarsN][t];
				}
				//---Vortex force
				AF_elem[m][k][0][t] = -cell_values[m][k][rho_id - 1][t] * cell_values[m][k][numVars + compVarsN + 13][t] * cell_values[m][k][numVars + compVarsN][t];
				AF_elem[m][k][1][t] = -cell_values[m][k][rho_id - 1][t] * cell_values[m][k][numVars + compVarsN + 14][t] * cell_values[m][k][numVars + compVarsN][t];
				//---Compressibility correction (Liu formulation - Volume integral term)
				AF_elem[m][k][20][t] = cell_values[m][k][numVars + 1][t] * cell_values[m][k][numVars+compVarsN+8][t] * cell_values[m][k][numVars + compVarsN][t];
				AF_elem[m][k][21][t] = cell_values[m][k][numVars + 1][t] * cell_values[m][k][numVars+compVarsN+9][t] * cell_values[m][k][numVars + compVarsN][t];
				//---Compressibility correction (original formulation)
				updateCompCorr(nodal_values,cell_values,rho_id,t,m,k,numVars,compVarsN,i_TOT,&kMax[0],NS,U_FORM,MOV_GRD,mrho_form,flux_scheme,near_ID,
								dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4,AF_elem);
				//---Compressibility correction (modified formulation)
				//------First term
				if (flags[m][k][0][t]==false) {
					updateFirstTerm(nodal_values,cell_values,rho_id,t,m,k,numVars,compVarsN,i_TOT,&kMax[0],NS,U_FORM,MOV_GRD,mrho_form,flux_scheme,near_ID,
									dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4,AF_elem);
				}
				if ((flags[m][k][0][t]==true) && (m>0) && (m<i_TOT-2) && (k<kMax[0]-2)) {	//Avoid borders as a check on neighbouring cells flag status will be performed.
																							//Elementary aerodynamic forces on cells at the border will anyhow never be used.
				//------Fifth term (modified for unsteady conditions)
					AF_elem[m][k][12][t]=(cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+13][t] -
										0.25*(aux_values_N[m][k][29][t]*nodal_values[m][k][vx_id-1][t]+aux_values_N[m+1][k][29][t]*nodal_values[m+1][k][vx_id-1][t]+
										aux_values_N[m+1][k+1][29][t]*nodal_values[m+1][k+1][vx_id-1][t]+aux_values_N[m][k+1][29][t]*nodal_values[m][k+1][vx_id-1][t])) *
										cell_values[m][k][numVars+compVarsN][t];
					AF_elem[m][k][13][t]=(cell_values[m][k][rho_id-1][t]*cell_values[m][k][numVars+compVarsN+14][t] -
										0.25*(aux_values_N[m][k][29][t]*nodal_values[m][k][vz_id-1][t]+aux_values_N[m+1][k][29][t]*nodal_values[m+1][k][vz_id-1][t]+
										aux_values_N[m+1][k+1][29][t]*nodal_values[m+1][k+1][vz_id-1][t]+aux_values_N[m][k+1][29][t]*nodal_values[m][k+1][vz_id-1][t])) *
										cell_values[m][k][numVars+compVarsN][t];
				//------Second term
					if (k != 0) {
						if (flags[m][k - 1][0][t] == false) {
							AF_elem[m][k][6][t] = AF_elem[m][k][6][t] - ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][vx_id - 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][vx_id - 1][t]) / 2) * (
									((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_1 +
									((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * (-dx_1));
						}
					}
					else {
						if ((flags[i_TOT-2-m][k][0][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
							AF_elem[m][k][6][t] = AF_elem[m][k][6][t] - ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][vx_id - 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][vx_id - 1][t]) / 2) * (
									((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_1 +
									((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * (-dx_1));
						} else if ((near_ID[m]) || (near_ID[m + 1])) {
							//Not effective in case of NS and fixed grid computation as V will be zero on the body surface
							AF_elem[m][k][6][t] = AF_elem[m][k][6][t] - ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][vx_id - 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][vx_id - 1][t]) / 2) * (
									((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_1 +
									((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * (-dx_1));
						}
					}
					if (flags[m+1][k][0][t]==false) {
						AF_elem[m][k][6][t]=AF_elem[m][k][6][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t]+
							nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vx_id-1][t])/2)*
							(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_2+
							((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_2));
					}
					if (flags[m][k+1][0][t]==false) {
						AF_elem[m][k][6][t]=AF_elem[m][k][6][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t]+
							nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vx_id-1][t])/2)*
							(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
							((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					}
					if (flags[m-1][k][0][t]==false) {
						AF_elem[m][k][6][t]=AF_elem[m][k][6][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vx_id-1][t]+
							nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vx_id-1][t])/2)*
							(((nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_4+
							((nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_4));
					}
					if (k != 0) {
						if (flags[m][k - 1][0][t] == false) {
							AF_elem[m][k][7][t] = AF_elem[m][k][7][t] - ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][vz_id - 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][vz_id - 1][t]) / 2) * (
									((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_1 +
									((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * (-dx_1));
						}
					}
					else {
						if ((flags[i_TOT - 2 - m][k][0][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
							AF_elem[m][k][7][t] = AF_elem[m][k][7][t] - ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][vz_id - 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][vz_id - 1][t]) / 2) * (
									((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_1 +
									((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * (-dx_1));
						} else if ((near_ID[m]) || (near_ID[m + 1])) {
							//Not effective in case of NS and fixed grid computation as V will be zero on the body surface
							AF_elem[m][k][7][t] = AF_elem[m][k][7][t] - ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][vz_id - 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][vz_id - 1][t]) / 2) * (
									((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_1 +
									((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * (-dx_1));
						}
					}
					if (flags[m+1][k][0][t]==false) {
						AF_elem[m][k][7][t]=AF_elem[m][k][7][t]-((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t]+
							nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][vz_id-1][t])/2)*
							(((nodal_values[m+1][k+1][vx_id-1][t]+nodal_values[m+1][k][vx_id-1][t])/2)*dz_2+
							((nodal_values[m+1][k+1][vz_id-1][t]+nodal_values[m+1][k][vz_id-1][t])/2)*(-dx_2));
					}
					if (flags[m][k+1][0][t]==false) {
						AF_elem[m][k][7][t]=AF_elem[m][k][7][t]-((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t]+
							nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][vz_id-1][t])/2)*
							(((nodal_values[m][k+1][vx_id-1][t]+nodal_values[m+1][k+1][vx_id-1][t])/2)*dz_3+
							((nodal_values[m][k+1][vz_id-1][t]+nodal_values[m+1][k+1][vz_id-1][t])/2)*(-dx_3));
					}
					if (flags[m-1][k][0][t]==false) {
						AF_elem[m][k][7][t]=AF_elem[m][k][7][t]-((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][vz_id-1][t]+
							nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][vz_id-1][t])/2)*
							(((nodal_values[m][k][vx_id-1][t]+nodal_values[m][k+1][vx_id-1][t])/2)*dz_4+
							((nodal_values[m][k][vz_id-1][t]+nodal_values[m][k+1][vz_id-1][t])/2)*(-dx_4));
					}
				//------Third term
					if (k != 0) {
						if (flags[m][k - 1][0][t] == false) {
							AF_elem[m][k][8][t] = AF_elem[m][k][8][t] + ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][numVars + 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][numVars + 1][t]) / 2) * dz_1;
						}
					}
					else {
						if ((flags[i_TOT - 2 - m][k][0][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
							AF_elem[m][k][8][t] = AF_elem[m][k][8][t] + ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][numVars + 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][numVars + 1][t]) / 2) * dz_1;
						} else if ((near_ID[m]) || (near_ID[m + 1])) {
							//Not effective in case of NS and fixed grid computation as V will be zero on the body surface
							AF_elem[m][k][8][t] = AF_elem[m][k][8][t] + ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][numVars + 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][numVars + 1][t]) / 2) * dz_1;
						}
					}
					if (flags[m+1][k][0][t]==false) {
						AF_elem[m][k][8][t]=AF_elem[m][k][8][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
							nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*dz_2;
					}
					if (flags[m][k+1][0][t]==false) {
						AF_elem[m][k][8][t]=AF_elem[m][k][8][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
							nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*dz_3;
					}
					if (flags[m-1][k][0][t]==false) {
						AF_elem[m][k][8][t]=AF_elem[m][k][8][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
							nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*dz_4;
					}
					if (k != 0) {
						if (flags[m][k - 1][0][t] == false) {
							AF_elem[m][k][9][t] = AF_elem[m][k][9][t] + ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][numVars + 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][numVars + 1][t]) / 2) * (-dx_1);
						}
					}
					else {
						if ((flags[i_TOT - 2 - m][k][0][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
							AF_elem[m][k][9][t] = AF_elem[m][k][9][t] + ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][numVars + 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][numVars + 1][t]) / 2) * (-dx_1);
						} else if ((near_ID[m]) || (near_ID[m + 1])) {
							//Not effective in case of NS and fixed grid computation as V will be zero on the body surface
							AF_elem[m][k][9][t] = AF_elem[m][k][9][t] + ((
								nodal_values[m + 1][k][rho_id - 1][t] * nodal_values[m + 1][k][numVars + 1][t] +
								nodal_values[m][k][rho_id - 1][t] * nodal_values[m][k][numVars + 1][t]) / 2) * (-dx_1);
						}
					}
					if (flags[m+1][k][0][t]==false) {
						AF_elem[m][k][9][t]=AF_elem[m][k][9][t]+((nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t]+
							nodal_values[m+1][k][rho_id-1][t]*nodal_values[m+1][k][numVars+1][t])/2)*(-dx_2);
					}
					if (flags[m][k+1][0][t]==false) {
						AF_elem[m][k][9][t]=AF_elem[m][k][9][t]+((nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t]+
							nodal_values[m+1][k+1][rho_id-1][t]*nodal_values[m+1][k+1][numVars+1][t])/2)*(-dx_3);
					}
					if (flags[m-1][k][0][t]==false) {
						AF_elem[m][k][9][t]=AF_elem[m][k][9][t]+((nodal_values[m][k][rho_id-1][t]*nodal_values[m][k][numVars+1][t]+
							nodal_values[m][k+1][rho_id-1][t]*nodal_values[m][k+1][numVars+1][t])/2)*(-dx_4);
					}
				//------Fourth term
					updateFourthTerm(nodal_values,aux_values_N,cell_values,t,m,k,numVars,compVarsN,i_TOT,flux_scheme,near_ID,
									rx_1,rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4,flags,AF_elem);
				}
				//---Parasite force (volume integrals formulation)
					if (!ParForm) {
						updateParasite(nodal_values,cell_values,aux_values,aux_values_N,flags,t,m,k,numVars,compVarsN,i_TOT,&kMax[0],NS,U_FORM,MOV_GRD,flux_scheme,smart_wds,near_ID,
										CdP_tsh,rho_inf,V_inf,chord,rx_1,rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4,AF_elem);
					} else {
						updateParasiteNew(nodal_values,cell_values,X_POLE,Z_POLE,rho_id,vx_id,vz_id,aux_values_N,flags,t,m,k,numVars,compVarsN,i_TOT,&kMax[0],NS,MOV_GRD,flux_scheme,smart_wds,near_ID,
										CdP_tsh,rho_inf,V_inf,chord,rx_1,rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4,AF_elem);
					}
					updateF_grad_rho(nodal_values,cell_values,aux_values,aux_values_N,flags,t,m,k,numVars,compVarsN,i_TOT,&kMax[0],NS,U_FORM,MOV_GRD,flux_scheme,smart_wds,near_ID,
									CdP_tsh,rho_inf,V_inf,rx_1,rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4,AF_elem);
				//FAR-FIELD elementary aerodynamic force using thermodynamic method
				if (thd_methods) {
					//DESTARAC & VAN DER VOOREN formula
					CdeP_thd_1= -V_inf * (
						(((nodal_values[m][k][vx_id - 1][t] + nodal_values[m + 1][k][vx_id - 1][t]) / 2) * dz_1 - ((nodal_values[m][k][vz_id - 1][t] + nodal_values[m + 1][k][vz_id - 1][t]) / 2) * dx_1) *
						((nodal_values[m][k][rho_id - 1][t] * aux_values_N[m][k][7][t] + nodal_values[m + 1][k][rho_id - 1][t] * aux_values_N[m + 1][k][7][t]) / 2) +
						(((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m + 1][k + 1][vx_id - 1][t]) / 2) * dz_2 - ((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m + 1][k + 1][vz_id - 1][t]) / 2) * dx_2) *
						((nodal_values[m + 1][k][rho_id - 1][t] * aux_values_N[m + 1][k][7][t] + nodal_values[m + 1][k + 1][rho_id - 1][t] * aux_values_N[m + 1][k + 1][7][t]) / 2) +
						(((nodal_values[m + 1][k + 1][vx_id - 1][t] + nodal_values[m][k + 1][vx_id - 1][t]) / 2) * dz_3 - ((nodal_values[m + 1][k + 1][vz_id - 1][t] + nodal_values[m][k + 1][vz_id - 1][t]) / 2) * dx_3) *
						((nodal_values[m + 1][k + 1][rho_id - 1][t] * aux_values_N[m + 1][k + 1][7][t] + nodal_values[m][k + 1][rho_id - 1][t] * aux_values_N[m][k + 1][7][t]) / 2) +
						(((nodal_values[m][k + 1][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_4 - ((nodal_values[m][k + 1][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * dx_4) *
						((nodal_values[m][k + 1][rho_id - 1][t] * aux_values_N[m][k + 1][7][t] + nodal_values[m][k][rho_id - 1][t] * aux_values_N[m][k][7][t]) / 2));	//+
						//dz_1*((aux_values_N[m][k][8] + aux_values_N[m + 1][k][8]) / 2) + dz_2 * ((aux_values_N[m+1][k][8] + aux_values_N[m + 1][k+1][8]) / 2)+
						//dz_3 * ((aux_values_N[m+1][k+1][8] + aux_values_N[m][k+1][8]) / 2)+ dz_4 * ((aux_values_N[m][k+1][8] + aux_values_N[m][k][8]) / 2));
					if (smart_wds) {
						if (((CdeP_thd_1> CdP_tsh* 0.5 * rho_inf * V_inf * V_inf * chord)|| (CdeP_thd_1 < -CdP_tsh*0.5*rho_inf*V_inf*V_inf*chord))&&flags[m][k][3][t]) {
							AF_elem[m][k][26][t] = CdeP_thd_1;
						} else if (flags[m][k][1][t]) {
							AF_elem[m][k][27][t] = CdeP_thd_1;
						} else {
							AF_elem[m][k][28][t] = CdeP_thd_1;
						}
					} else {
						if (flags[m][k][3][t]) {
							AF_elem[m][k][26][t] = CdeP_thd_1;
						} else if (flags[m][k][1][t]) {
							AF_elem[m][k][27][t] = CdeP_thd_1;
						} else {
							AF_elem[m][k][28][t] = CdeP_thd_1;
						}
					}
					//PAPARONE-TOGNACCINI formula (Entropy drag)
					CdeP_thd_2= -V_inf * (
						(((nodal_values[m][k][vx_id - 1][t] + nodal_values[m + 1][k][vx_id - 1][t]) / 2) * dz_1 - ((nodal_values[m][k][vz_id - 1][t] + nodal_values[m + 1][k][vz_id - 1][t]) / 2) * dx_1) *
						((nodal_values[m][k][rho_id - 1][t] * aux_values_N[m][k][6][t] + nodal_values[m + 1][k][rho_id - 1][t] * aux_values_N[m + 1][k][6][t]) / 2) +
						(((nodal_values[m + 1][k][vx_id - 1][t] + nodal_values[m + 1][k + 1][vx_id - 1][t]) / 2) * dz_2 - ((nodal_values[m + 1][k][vz_id - 1][t] + nodal_values[m + 1][k + 1][vz_id - 1][t]) / 2) * dx_2) *
						((nodal_values[m + 1][k][rho_id - 1][t] * aux_values_N[m + 1][k][6][t] + nodal_values[m + 1][k + 1][rho_id - 1][t] * aux_values_N[m + 1][k + 1][6][t]) / 2) +
						(((nodal_values[m + 1][k + 1][vx_id - 1][t] + nodal_values[m][k + 1][vx_id - 1][t]) / 2) * dz_3 - ((nodal_values[m + 1][k + 1][vz_id - 1][t] + nodal_values[m][k + 1][vz_id - 1][t]) / 2) * dx_3) *
						((nodal_values[m + 1][k + 1][rho_id - 1][t] * aux_values_N[m + 1][k + 1][6][t] + nodal_values[m][k + 1][rho_id - 1][t] * aux_values_N[m][k + 1][6][t]) / 2) +
						(((nodal_values[m][k + 1][vx_id - 1][t] + nodal_values[m][k][vx_id - 1][t]) / 2) * dz_4 - ((nodal_values[m][k + 1][vz_id - 1][t] + nodal_values[m][k][vz_id - 1][t]) / 2) * dx_4) *
						((nodal_values[m][k + 1][rho_id - 1][t] * aux_values_N[m][k + 1][6][t] + nodal_values[m][k][rho_id - 1][t] * aux_values_N[m][k][6][t]) / 2));
					if (smart_wds) {
						if (((CdeP_thd_2 > CdP_tsh * 0.5 * rho_inf * V_inf * V_inf * chord) || (CdeP_thd_2 < -CdP_tsh * 0.5 * rho_inf * V_inf * V_inf * chord)) && flags[m][k][3][t]) {
							AF_elem[m][k][29][t] = CdeP_thd_2;
						} else if (flags[m][k][1][t]) {
							AF_elem[m][k][30][t] = CdeP_thd_2;
						} else {
							AF_elem[m][k][31][t] = CdeP_thd_2;
							flags[m][k][3][t] = false;
						}
					} else {
						if (flags[m][k][3][t]) {
							AF_elem[m][k][29][t] = CdeP_thd_2;
						} else if (flags[m][k][1][t]) {
							AF_elem[m][k][30][t] = CdeP_thd_2;
						} else {
							AF_elem[m][k][31][t] = CdeP_thd_2;
						}
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of elementary aerodynamic forces calculation.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Compute spatial gradients of (rho*dV/dt) or (rho*dV/dt') at grid nodes.
	start = chrono::steady_clock::now();
	if ((grad_scheme==1)&&(!FT_FORM)) {	//Only compute in case they are needed (FT_FORM==false) and the FD gradient scheme is selected (grad_scheme==1)
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k <kMax[0]; k++) {
				for (int m = 0; m <i_TOT; m++) {
					//Local Jacobian
					double J11=0, J12=0, J21=0, J22=0, ddm=0, ddk=0, drhoVtz_dx=0, drhoVtx_dz=0;
					if (m==0) {
						J11=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m+1][k][0][t]-0.5*nodal_values[m+2][k][0][t];
						J12=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m+1][k][2][t]-0.5*nodal_values[m+2][k][2][t];
					} else if (m==i_TOT-1) {
						J11=0.5*nodal_values[m-2][k][0][t]-2.0*nodal_values[m-1][k][0][t]+1.5*nodal_values[m][k][0][t];
						J12=0.5*nodal_values[m-2][k][2][t]-2.0*nodal_values[m-1][k][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J11=0.5*nodal_values[m+1][k][0][t]-0.5*nodal_values[m-1][k][0][t];
						J12=0.5*nodal_values[m+1][k][2][t]-0.5*nodal_values[m-1][k][2][t];
					}
					if (k==0) {
						J21=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k+2][0][t];
						J22=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k+2][2][t];
					} else if (k==kMax[0]-1) {
						J21=0.5*nodal_values[m][k-2][0][t]-2.0*nodal_values[m][k-1][0][t]+1.5*nodal_values[m][k][0][t];
						J22=0.5*nodal_values[m][k-2][2][t]-2.0*nodal_values[m][k-1][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J21=0.5*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k-1][0][t];
						J22=0.5*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k-1][2][t];
					}
					double detJ=J11*J22-J12*J21;
					if (U_FORM==0) {				//Index gradients of (rho*dV/dt)
						if (m==0) {
							ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][27][t]+
								 2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][27][t]
								-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][27][t];
						} else if (m==i_TOT-1) {
							ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][27][t]
								-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][27][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][27][t];
						} else {
							ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][27][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][27][t];
						}
						if (k==0) {
							ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][27][t]
								+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][27][t]
								-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][27][t];
						} else if (k==kMax[0]-1) {
							ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][27][t]
								-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][27][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][27][t];
						} else {
							ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][27][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][27][t];
						}
						//Spatial gradients
						drhoVtx_dz = (-J21*ddm + J11*ddk) / detJ;
						if (m==0) {
							ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][28][t]+
								 2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][28][t]
								-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][28][t];
						} else if (m==i_TOT-1) {
							ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][28][t]
								-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][28][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][28][t];
						} else {
							ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][28][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][28][t];
						}
						if (k==0) {
							ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][28][t]
								+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][28][t]
								-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][28][t];
						} else if (k==kMax[0]-1) {
							ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][28][t]
								-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][28][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][28][t];
						} else {
							ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][28][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][28][t];
						}
						//Spatial gradients
						drhoVtz_dx = (J22*ddm - J12*ddk) / detJ;
						//TEMPORARILY store curl values in the nodal position of Y coordinates (never used)
						nodal_values[m][k][1][t] = drhoVtx_dz - drhoVtz_dx;
					} else if (U_FORM==1) {			//Index gradients of (rho*dV'/dt')
						if (m==0) {
							ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t]+
								 2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t]
								-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][25][t];
						} else if (m==i_TOT-1) {
							ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][25][t]
								-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][25][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t];
						} else {
							ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][25][t];
						}
						if (k==0) {
							ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t]
								+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t]
								-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][25][t];
						} else if (k==kMax[0]-1) {
							ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][25][t]
								-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][25][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t];
						} else {
							ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][25][t];
						}
						//Spatial gradients
						drhoVtx_dz = (-J21*ddm + J11*ddk) / detJ;
						if (m==0) {
							ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t]+
								 2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t]
								-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][26][t];
						} else if (m==i_TOT-1) {
							ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][26][t]
								-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][26][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t];
						} else {
							ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][26][t];
						}
						if (k==0) {
							ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t]
								+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t]
								-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][26][t];
						} else if (k==kMax[0]-1) {
							ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][26][t]
								-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][26][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t];
						} else {
							ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][26][t];
						}
						//Spatial gradients
						drhoVtz_dx = (J22*ddm - J12*ddk) / detJ;
						//TEMPORARILY store curl values in the nodal position of Y coordinates (never used)
						nodal_values[m][k][1][t] = drhoVtx_dz - drhoVtz_dx;
					} else if (U_FORM==2) {			//Index gradients of (rho*dV/dt')
						if (m==0) {
							ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t]+
								 2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t]
								-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][25][t];
						} else if (m==i_TOT-1) {
							ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][25][t]
								-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][25][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t];
						} else {
							ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][25][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][25][t];
						}
						if (k==0) {
							ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t]
								+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t]
								-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][25][t];
						} else if (k==kMax[0]-1) {
							ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][25][t]
								-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][25][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][25][t];
						} else {
							ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][25][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][25][t];
						}
						//Spatial gradients
						drhoVtx_dz = (-J21*ddm + J11*ddk) / detJ;
						if (m==0) {
							ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t]+
								 2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t]
								-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][26][t];
						} else if (m==i_TOT-1) {
							ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][26][t]
								-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][26][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t];
						} else {
							ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][26][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][26][t];
						}
						if (k==0) {
							ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t]
								+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t]
								-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][26][t];
						} else if (k==kMax[0]-1) {
							ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][26][t]
								-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][26][t]
								+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][26][t];
						} else {
							ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][26][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][26][t];
						}
						//Spatial gradients
						drhoVtz_dx = (J22*ddm - J12*ddk) / detJ;
						//TEMPORARILY store curl values in the nodal position of Y coordinates (never used)
						nodal_values[m][k][1][t] = drhoVtx_dz - drhoVtz_dx;
					}			
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of computing spatial gradients of of (rho*dV/dt) or (rho*dV/dt') at grid nodes.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Compute spatial gradients of (rho*grad(kt)) at grid nodes.
	start = chrono::steady_clock::now();
	if ((grad_scheme==1)&&(Lamb_form>0)) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k <kMax[0]; k++) {
				for (int m = 0; m <i_TOT; m++) {
					//Local Jacobian
					double J11=0, J12=0, J21=0, J22=0, ddm=0, ddk=0, drhoGradXkt_dz=0, drhoGradZkt_dx=0;
					if (m==0) {
						J11=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m+1][k][0][t]-0.5*nodal_values[m+2][k][0][t];
						J12=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m+1][k][2][t]-0.5*nodal_values[m+2][k][2][t];
					} else if (m==i_TOT-1) {
						J11=0.5*nodal_values[m-2][k][0][t]-2.0*nodal_values[m-1][k][0][t]+1.5*nodal_values[m][k][0][t];
						J12=0.5*nodal_values[m-2][k][2][t]-2.0*nodal_values[m-1][k][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J11=0.5*nodal_values[m+1][k][0][t]-0.5*nodal_values[m-1][k][0][t];
						J12=0.5*nodal_values[m+1][k][2][t]-0.5*nodal_values[m-1][k][2][t];
					}
					if (k==0) {
						J21=-1.5*nodal_values[m][k][0][t]+2.0*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k+2][0][t];
						J22=-1.5*nodal_values[m][k][2][t]+2.0*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k+2][2][t];
					} else if (k==kMax[0]-1) {
						J21=0.5*nodal_values[m][k-2][0][t]-2.0*nodal_values[m][k-1][0][t]+1.5*nodal_values[m][k][0][t];
						J22=0.5*nodal_values[m][k-2][2][t]-2.0*nodal_values[m][k-1][2][t]+1.5*nodal_values[m][k][2][t];
					} else {
						J21=0.5*nodal_values[m][k+1][0][t]-0.5*nodal_values[m][k-1][0][t];
						J22=0.5*nodal_values[m][k+1][2][t]-0.5*nodal_values[m][k-1][2][t];
					}
					double detJ=J11*J22-J12*J21;
					//Indicial gradients
					if (m==0) {
						ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][30][t]+
								2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][30][t]
							-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][30][t];
					} else if (m==i_TOT-1) {
						ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][30][t]
							-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][30][t]
							+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][30][t];
					} else {
						ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][30][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][30][t];
					}
					if (k==0) {
						ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][30][t]
							+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][30][t]
							-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][30][t];
					} else if (k==kMax[0]-1) {
						ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][30][t]
							-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][30][t]
							+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][30][t];
					} else {
						ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][30][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][30][t];
					}
					//Spatial gradients
					drhoGradXkt_dz = (-J21*ddm + J11*ddk) / detJ;
					//Indicial gradients
					if (m==0) {
						ddm=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][31][t]+
								2.0*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][31][t]
							-0.5*nodal_values[m+2][k][rho_id-1][t]*aux_values_N[m+2][k][31][t];
					} else if (m==i_TOT-1) {
						ddm= 0.5*nodal_values[m-2][k][rho_id-1][t]*aux_values_N[m-2][k][31][t]
							-2.0*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][31][t]
							+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][31][t];
					} else {
						ddm=0.5*nodal_values[m+1][k][rho_id-1][t]*aux_values_N[m+1][k][31][t]-0.5*nodal_values[m-1][k][rho_id-1][t]*aux_values_N[m-1][k][31][t];
					}
					if (k==0) {
						ddk=-1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][31][t]
							+2.0*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][31][t]
							-0.5*nodal_values[m][k+2][rho_id-1][t]*aux_values_N[m][k+2][31][t];
					} else if (k==kMax[0]-1) {
						ddk= 0.5*nodal_values[m][k-2][rho_id-1][t]*aux_values_N[m][k-2][31][t]
							-2.0*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][31][t]
							+1.5*nodal_values[m][k][rho_id-1][t] * aux_values_N[m][k][31][t];
					} else {
						ddk=0.5*nodal_values[m][k+1][rho_id-1][t]*aux_values_N[m][k+1][31][t]-0.5*nodal_values[m][k-1][rho_id-1][t]*aux_values_N[m][k-1][31][t];
					}
					//Spatial gradients
					drhoGradZkt_dx = (J22*ddm - J12*ddk) / detJ;
					//Store curl values in the auxiliary nodal position [32] (which was for TEMPORARY use, storing total enthalpy gradient along X on grid nodes)
					aux_values_N[m][k][32][t] = drhoGradXkt_dz - drhoGradZkt_dx;
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of computing spatial gradients of (rho*grad(kt)) at grid nodes..   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Update volume integrand of unsteady term, in case of FD gradient scheme selected
	start = chrono::steady_clock::now();
	if ((grad_scheme==1)&&(!FT_FORM)) {
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k <kMax[0]-1; k++) {
				for (int m = 0; m <i_TOT-1; m++) {
					AF_elem[m][k][35][t] = cell_values[m][k][numVars + compVarsN + 16][t] * 0.25 *
						(nodal_values[m][k][1][t]+nodal_values[m+1][k][1][t]+nodal_values[m+1][k+1][1][t]+nodal_values[m][k+1][1][t]) *
						cell_values[m][k][numVars + compVarsN][t];
					AF_elem[m][k][36][t] = -cell_values[m][k][numVars + compVarsN + 15][t] * 0.25 *
						(nodal_values[m][k][1][t]+nodal_values[m+1][k][1][t]+nodal_values[m+1][k+1][1][t]+nodal_values[m][k+1][1][t]) *
						cell_values[m][k][numVars + compVarsN][t];
					if (Lamb_form > 0) {	//Add turbulent term
						AF_elem[m][k][35][t] = AF_elem[m][k][35][t] + cell_values[m][k][numVars + compVarsN + 16][t] *
												(0.25*(aux_values_N[m][k][32][t]+aux_values_N[m+1][k][32][t]+aux_values_N[m+1][k+1][32][t]+aux_values_N[m][k+1][32][t])) *
												cell_values[m][k][numVars + compVarsN][t];
						AF_elem[m][k][36][t] = AF_elem[m][k][36][t] -cell_values[m][k][numVars + compVarsN + 15][t] *
												(0.25*(aux_values_N[m][k][32][t]+aux_values_N[m+1][k][32][t]+aux_values_N[m+1][k+1][32][t]+aux_values_N[m][k+1][32][t])) *
												cell_values[m][k][numVars + compVarsN][t];
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of updating volume integrand of unsteady term.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//---Conditionally restore Y coordinate values.
	start = chrono::steady_clock::now();
	if (grad_scheme==1) {	//Only restore in cases they were altered (search for "TEMPORARILY store..." comments)
		for (int t=0; t<Nt; t++) {
			for (int k = 0; k <kMax[0]; k++) {
				for (int m = 0; m <i_TOT; m++) {
					nodal_values[m][k][1][t]=0;
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of restoring Y coordinate values.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Smart flagging in the viscous wake region
	start = chrono::steady_clock::now();
	if (dom_decomp) {
		/* for (int t=0; t<Nt; t++) {
			for (int m = m_te-10; m <i_TOT-2-m_te+1+10; m++) {
				int k=112;
				//double Cdp_vrt=(AF_elem[m][k][14][t] + AF_elem[m][k][15][t] + AF_elem[m][k][17][t] + AF_elem[m][k][18][t] + AF_elem[m][k][33][t] + AF_elem[m][k][34][t] +
				//				AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t])/ (0.5 * rho_inf * pow(V_inf, 2) * chord);
				double Cdp_vrt=(AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t])/ (0.5 * rho_inf * pow(V_inf, 2) * chord);
				while ((k>0)&&(abs(Cdp_vrt)<CdP_tsh)) {
					k--;
					Cdp_vrt=(AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t])/ (0.5 * rho_inf * pow(V_inf, 2) * chord);
				}
				for (int kk=0;kk<=k;kk++) {
					flags[m][kk][1][t]=1;
					flags[m][kk][2][t]=0;
					flags[m][kk][3][t]=0;
				}
			}
		} */
		/* for (int t=0; t<Nt; t++) {
			for (int m = 0; m <m_te; m++) {
				int k=kMax[0]-2-10;
				double Cdp_vrt=(AF_elem[m][k][14][t] + AF_elem[m][k][15][t] + AF_elem[m][k][17][t] + AF_elem[m][k][18][t] + AF_elem[m][k][33][t] + AF_elem[m][k][34][t] +
								AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t])/ (0.5 * rho_inf * pow(V_inf, 2) * chord);
				while ((k>0)&&(abs(Cdp_vrt)<CdP_tsh)) {
					k--;
					Cdp_vrt=(AF_elem[m][k][14][t] + AF_elem[m][k][15][t] + AF_elem[m][k][17][t] + AF_elem[m][k][18][t] + AF_elem[m][k][33][t] + AF_elem[m][k][34][t] +
							AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t])/ (0.5 * rho_inf * pow(V_inf, 2) * chord);
				}
				for (int kk=0;kk<=k;kk++) {
					flags[m][kk][1][t]=1;
					flags[m][kk][2][t]=0;
					flags[m][kk][3][t]=0;
				}
			}
			for (int m = i_TOT-2-m_te+1; m <=i_TOT-2; m++) {
				int k=kMax[0]-2-10;
				double Cdp_vrt=(AF_elem[m][k][14][t] + AF_elem[m][k][15][t] + AF_elem[m][k][17][t] + AF_elem[m][k][18][t] + AF_elem[m][k][33][t] + AF_elem[m][k][34][t] +
								AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t])/ (0.5 * rho_inf * pow(V_inf, 2) * chord);
				while ((k>0)&&(abs(Cdp_vrt)<CdP_tsh)) {
					k--;
					Cdp_vrt=(AF_elem[m][k][14][t] + AF_elem[m][k][15][t] + AF_elem[m][k][17][t] + AF_elem[m][k][18][t] + AF_elem[m][k][33][t] + AF_elem[m][k][34][t] +
							AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t])/ (0.5 * rho_inf * pow(V_inf, 2) * chord);
				}
				//cout << "k=" << k << "\n";
				for (int kk=0;kk<=k;kk++) {
					flags[m][kk][1][t]=1;
					flags[m][kk][2][t]=0;
					flags[m][kk][3][t]=0;
				}
			}
		} */
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of smart flagging in the viscous wake region.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Surface integrals on wave region boundary
	start = chrono::steady_clock::now();
	if (dom_decomp) {
		//
		//WARNING:   the Lamb vector circulation moment and the ONERA compressibility correction are BOTH integrated using the variable AF_elem[m][k][32][t]
		//
		/* for (int t=0; t<Nt; t++) {
			for (int k = 0; k < kMax[0] - 2; k++) {
				for (int m = 1; m < i_TOT - 2; m++) {			
					if (flags[m][k][2][t]) {
						cell_metrics(nodal_values, cell_values, NS, U_FORM, MOV_GRD, t, m, k, numVars, compVarsN, near_ID, i_TOT, &kMax[0], grad_scheme, flux_scheme, ParForm, X_POLE, Z_POLE,
										rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS);
						if (flux_scheme) {
							if (k!=0) {
								if (flags[m][k - 1][2][t] == false) {
									//Lamb vector circulation moment
		 							AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
										rz_1 * ((nodal_values[m][k][rho_id - 1][t] + nodal_values[m + 1][k][rho_id - 1][t]) / 2) * (
										dz_1 * (a_1 * cell_values[m][k-1][numVars+compVarsN+14][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+14][t]) +
										dx_1 * (a_1 * cell_values[m][k-1][numVars+compVarsN+13][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+13][t]));
									//ONERA compressibility correction term
									AF_elem[m][k][32][t] = AF_elem[m][k][32][t] + 
													0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
													rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
													(dz_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_1 * (a_1 * cell_values[m][k - 1][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
								}
							} else {
								if ((flags[i_TOT - 2 - m][k][2][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
									//Lamb vector circulation moment
		 							AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
										rz_1 * ((nodal_values[i_TOT - 2 - m][k][rho_id - 1][t] + nodal_values[m + 1][k][rho_id - 1][t]) / 2) * (
										dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k-1][numVars+compVarsN+14][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+14][t]) +
										dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k-1][numVars+compVarsN+13][t] + (1 - a_1) * cell_values[m][k][numVars+compVarsN+13][t]));
									//ONERA compressibility correction term
									AF_elem[m][k][32][t] = AF_elem[m][k][32][t] + 
													0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
													rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m + 1][k][numVars + 1][t]) / 2) *
													(dz_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 9][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 9][t]) +
													dx_1 * (a_1 * cell_values[i_TOT - 2 - m][k][numVars + compVarsN + 8][t] + (1 - a_1) * cell_values[m][k][numVars + compVarsN + 8][t]));
								}
							}
							if ((flags[m + 1][k][2][t] == false)||(m==i_TOT-3)) {
								//Lamb vector circulation moment
		 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
									rz_2 * ((nodal_values[m + 1][k][rho_id - 1][t] + nodal_values[m + 1][k + 1][rho_id - 1][t]) / 2) * (
									dz_2 * (a_2 * cell_values[m+1][k][numVars+compVarsN+14][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+14][t]) +
									dx_2 * (a_2 * cell_values[m+1][k][numVars+compVarsN+13][t] + (1 - a_2) * cell_values[m][k][numVars+compVarsN+13][t]));
								//ONERA compressibility correction term
								AF_elem[m][k][32][t] = AF_elem[m][k][32][t] + 
										0.5*pow(V_inf,2)*(((nodal_values[m + 1][k][rho_id-1][t]+nodal_values[m + 1][k + 1][rho_id-1][t])/2)-rho_inf)*dz_2 +
										rz_2 * ((nodal_values[m + 1][k][numVars + 1][t] + nodal_values[m + 1][k + 1][numVars + 1][t]) / 2) *
										(dz_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 9][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 9][t]) +
										dx_2 * (a_2 * cell_values[m + 1][k][numVars + compVarsN + 8][t] + (1 - a_2) * cell_values[m][k][numVars + compVarsN + 8][t]));
							}
							if ((flags[m][k + 1][2][t] == false)||(k==kMax[0]-3)) {
								//Lamb vector circulation moment
		 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
									rz_3 * ((nodal_values[m + 1][k + 1][rho_id - 1][t] + nodal_values[m][k + 1][rho_id - 1][t]) / 2) * (
									dz_3 * (a_3 * cell_values[m][k+1][numVars+compVarsN+14][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+14][t]) +
									dx_3 * (a_3 * cell_values[m][k+1][numVars+compVarsN+13][t] + (1 - a_3) * cell_values[m][k][numVars+compVarsN+13][t]));
								//ONERA compressibility correction term
								AF_elem[m][k][32][t] = AF_elem[m][k][32][t] + 
										0.5*pow(V_inf,2)*(((nodal_values[m + 1][k + 1][rho_id-1][t]+nodal_values[m][k + 1][rho_id-1][t])/2)-rho_inf)*dz_3 +
										rz_3 * ((nodal_values[m + 1][k + 1][numVars + 1][t] + nodal_values[m][k + 1][numVars + 1][t]) / 2) *
										(dz_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 9][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 9][t]) +
										dx_3 * (a_3 * cell_values[m][k + 1][numVars + compVarsN + 8][t] + (1 - a_3) * cell_values[m][k][numVars + compVarsN + 8][t]));
							}
							if ((flags[m - 1][k][2][t] == false)||(m==1)) {
								//Lamb vector circulation moment
		 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
									rz_4 * ((nodal_values[m][k + 1][rho_id - 1][t] + nodal_values[m][k][rho_id - 1][t]) / 2) * (
									dz_4 * (a_4 * cell_values[m-1][k][numVars+compVarsN+14][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+14][t]) +
									dx_4 * (a_4 * cell_values[m-1][k][numVars+compVarsN+13][t] + (1 - a_4) * cell_values[m][k][numVars+compVarsN+13][t]));
								//ONERA compressibility correction term
								AF_elem[m][k][32][t] = AF_elem[m][k][32][t] + 
										0.5*pow(V_inf,2)*(((nodal_values[m][k + 1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
										rz_4 * ((nodal_values[m][k + 1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
										(dz_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 9][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 9][t]) +
										dx_4 * (a_4 * cell_values[m - 1][k][numVars + compVarsN + 8][t] + (1 - a_4) * cell_values[m][k][numVars + compVarsN + 8][t]));
							}
						} else {
							if (k!=0) {
								if (flags[m][k - 1][2][t] == false) {
									//Lamb vector circulation moment
			 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
										rz_1 * ((nodal_values[m][k][rho_id - 1][t] + nodal_values[m + 1][k][rho_id - 1][t]) / 2) * (
										dz_1 * ((aux_values_N[m][k][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
										dx_1 * ((aux_values_N[m][k][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
									//ONERA compressibility correction term
									AF_elem[m][k][32][t] = AF_elem[m][k][32][t] +
											0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
											rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) * (
											dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
											dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
								}
							} else {
								if ((flags[i_TOT - 2 - m][k][2][t] == false) && (!near_ID[m]) && (!near_ID[m+1])) {
									//Lamb vector circulation moment
			 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
										rz_1 * ((nodal_values[m][k][rho_id - 1][t] + nodal_values[m + 1][k][rho_id - 1][t]) / 2) * (
										dz_1 * ((aux_values_N[m][k][20][t] + aux_values_N[m + 1][k][20][t]) / 2) +
										dx_1 * ((aux_values_N[m][k][19][t] + aux_values_N[m + 1][k][19][t]) / 2));
									//ONERA compressibility correction term
									AF_elem[m][k][32][t] = AF_elem[m][k][32][t] +
											0.5*pow(V_inf,2)*(((nodal_values[m][k][rho_id-1][t]+nodal_values[m+1][k][rho_id-1][t])/2)-rho_inf)*dz_1 +
											rz_1 * ((nodal_values[m][k][numVars + 1][t] + nodal_values[m+1][k][numVars + 1][t]) / 2) * (
											dz_1 * ((aux_values_N[m][k][18][t] + aux_values_N[m+1][k][18][t]) / 2) +
											dx_1 * ((aux_values_N[m][k][17][t] + aux_values_N[m+1][k][17][t]) / 2));
								}
							}
							if ((flags[m + 1][k][2][t] == false)||(m==i_TOT-3)) {
								//Lamb vector circulation moment
		 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
									rz_2 * ((nodal_values[m + 1][k][rho_id - 1][t] + nodal_values[m + 1][k + 1][rho_id - 1][t]) / 2) * (
									dz_2 * ((aux_values_N[m + 1][k][20][t] + aux_values_N[m + 1][k + 1][20][t]) / 2) +
									dx_2 * ((aux_values_N[m + 1][k][19][t] + aux_values_N[m + 1][k + 1][19][t]) / 2));
								//ONERA compressibility correction term
								AF_elem[m][k][32][t] = AF_elem[m][k][32][t] +
										0.5*pow(V_inf,2)*(((nodal_values[m+1][k][rho_id-1][t]+nodal_values[m+1][k+1][rho_id-1][t])/2)-rho_inf)*dz_2 +
										rz_2 * ((nodal_values[m+1][k][numVars + 1][t] + nodal_values[m+1][k+1][numVars + 1][t]) / 2) *
										(dz_2 * ((aux_values_N[m+1][k][18][t] + aux_values_N[m+1][k+1][18][t]) / 2) +
										dx_2 * ((aux_values_N[m+1][k][17][t] + aux_values_N[m+1][k+1][17][t]) / 2));
							}
							if ((flags[m][k + 1][2][t] == false)||(k==kMax[0]-3)) {
								//Lamb vector circulation moment
		 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
									rz_3 * ((nodal_values[m + 1][k + 1][rho_id - 1][t] + nodal_values[m][k + 1][rho_id - 1][t]) / 2) * (
									dz_3 * ((aux_values_N[m + 1][k + 1][20][t] + aux_values_N[m][k + 1][20][t]) / 2) +
									dx_3 * ((aux_values_N[m + 1][k + 1][19][t] + aux_values_N[m][k + 1][19][t]) / 2));
								//ONERA compressibility correction term
								AF_elem[m][k][32][t] = AF_elem[m][k][32][t] +
										0.5*pow(V_inf,2)*(((nodal_values[m+1][k+1][rho_id-1][t]+nodal_values[m][k+1][rho_id-1][t])/2)-rho_inf)*dz_3 +
										rz_3 * ((nodal_values[m+1][k+1][numVars + 1][t] + nodal_values[m][k+1][numVars + 1][t]) / 2) *
										(dz_3 * ((aux_values_N[m+1][k+1][18][t] + aux_values_N[m][k+1][18][t]) / 2) +
										dx_3 * ((aux_values_N[m+1][k+1][17][t] + aux_values_N[m][k+1][17][t]) / 2));
							}
							if ((flags[m - 1][k][2][t] == false)||(m==1)) {
								//Lamb vector circulation moment
		 						AF_elem[m][k][32][t] = AF_elem[m][k][32][t] -
									rz_4 * ((nodal_values[m][k + 1][rho_id - 1][t] + nodal_values[m][k][rho_id - 1][t]) / 2) * (
									dz_4 * ((aux_values_N[m][k + 1][20][t] + aux_values_N[m][k][20][t]) / 2) +
									dx_4 * ((aux_values_N[m][k + 1][19][t] + aux_values_N[m][k][19][t]) / 2));
								//ONERA compressibility correction term
								AF_elem[m][k][32][t] = AF_elem[m][k][32][t] +
										0.5*pow(V_inf,2)*(((nodal_values[m][k+1][rho_id-1][t]+nodal_values[m][k][rho_id-1][t])/2)-rho_inf)*dz_4 +
										rz_4 * ((nodal_values[m][k+1][numVars + 1][t] + nodal_values[m][k][numVars + 1][t]) / 2) *
										(dz_4 * ((aux_values_N[m][k+1][18][t] + aux_values_N[m][k][18][t]) / 2) +
										dx_4 * ((aux_values_N[m][k+1][17][t] + aux_values_N[m][k][17][t]) / 2));
							}
						}
					}
			// 		//if (flags[m][k][3]) {
			// 		//	CdP_vrt_thd = CdP_vrt_thd - (((aux_values_N[m][k][0] + aux_values_N[m + 1][k][0]) / 2) * dz_1 +
			// 		//		((aux_values_N[m + 1][k][0] + aux_values_N[m + 1][k + 1][0]) / 2) * dz_2 +
			// 		//		((aux_values_N[m][k + 1][0] + aux_values_N[m + 1][k + 1][0]) / 2) * dz_3 +
			// 		//		((aux_values_N[m][k][0] + aux_values_N[m][k + 1][0]) / 2) * dz_4);
			// 		//	CdP_vrt_thd = CdP_vrt_thd + (
			// 		//		(rx_1 * dz_1 - rz_1 * dx_1) * ((aux_values_N[m][k][1] + aux_values_N[m + 1][k][1]) / 2) +
			// 		//		(rx_2 * dz_2 - rz_2 * dx_2) * ((aux_values_N[m + 1][k][1] + aux_values_N[m + 1][k + 1][1]) / 2) +
			// 		//		(rx_3 * dz_3 - rz_3 * dx_3) * ((aux_values_N[m][k + 1][1] + aux_values_N[m + 1][k + 1][1]) / 2) +
			// 		//		(rx_4 * dz_4 - rz_4 * dx_4) * ((aux_values_N[m][k][1] + aux_values_N[m][k + 1][1]) / 2));
			// 		//}
 				}
			}
		} */
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of surface integrals calculation on wave region boundary.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Aerodynamic force non-dimensionalization
	start = chrono::steady_clock::now();
	for (int t=0; t<Nt; t++) {
		Clf_nf[t] = Clf_nf[t] / (0.5 * rho_inf * pow(V_inf, 2) * chord);
		Clp_nf[t] = Clp_nf[t] / (0.5 * rho_inf * pow(V_inf, 2) * chord);
		Cdf_nf[t] = Cdf_nf[t] / (0.5 * rho_inf * pow(V_inf, 2) * chord);
		Cdp_nf[t] = Cdp_nf[t] / (0.5 * rho_inf * pow(V_inf, 2) * chord);
		for (int j = 0; j <numAF; j++) {
			for (int k = 0; k <kMax[0]-1; k++) {
				for (int m = 0; m <i_TOT-1; m++) {
					if (AF_elem[m][k][j][t]!=0) {
						AF_elem[m][k][j][t]=AF_elem[m][k][j][t]/(0.5*rho_inf*pow(V_inf,2)*chord);
					}
				}
			}
		}
	}
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of aerodynamic force non-dimensionalization.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Compute aerodynamic force components
	vector < vector < double > > AERO(numAF+numAF_ADD, vector < double > (Nt, 0));
	vector < vector < vector < vector < double > > > > AF;
	if ((k_SEL>0)&&(m_SEL>0)) {
		//Compute aerodynamic force components for the selected integration domain
		start = chrono::steady_clock::now();
		for (int t=0; t<Nt; t++) {
			compute_AF(nodal_values, cell_values, aux_values, aux_values_N, flags, AF_elem, t, int(m_SEL), int(k_SEL), numVars, compVarsN, numAF, numAF_ADD, i_TOT, &kMax[0], NS, flux_scheme,
					grad_scheme, ParForm, limit_integr, near_ID, rho_inf, V_inf, chord, rho_id, vx_id, vz_id, X_POLE, Z_POLE, U_FORM, FT_FORM, x_CoR, z_CoR, ampl_x, freq_x, t0_x, ampl_z, freq_z, t0_z,
					ampl_a, freq_a, t0_a, MOV_GRD, time, rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4, dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS, AERO);
		}
		if (debug_mode) {
			end = chrono::steady_clock::now();
			cout << "END of computing aerodynamic force components for the selected integration domain.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
		}
	} else {
		AERO.clear();
		AF = vector < vector < vector < vector < double > > > >(floor((i_TOT-1)/2)+1, vector < vector < vector < double > > > (kMax[0], vector < vector < double > > (numAF+numAF_ADD, vector < double > (Nt,0))));
		//Include near-field aerodynamic force values in the AF data structure
		start = chrono::steady_clock::now();
		for (int t=0; t<Nt; t++) {
			int mD=0;
			for (int kD = 0; kD <kMax[0]; kD++) {
				AF[mD][kD][0][t]=Clp_nf[t]+Clf_nf[t];        //Set 1st aerodynamic force component at near field lift coeff. value at grid far-field boundary
				AF[mD][kD][1][t]=Cdp_nf[t]+Cdf_nf[t];        //Set 2nd aerodynamic force component at near field drag coeff. value at grid far-field boundary
				AF[mD][kD][2][t]=Clf_nf[t];        			//Set 3rd aerodynamic force component at near field friction lift coeff. value at grid far-field boundary
				AF[mD][kD][3][t]=Cdf_nf[t];        			//Set 4th aerodynamic force component at near field friction drag coeff. value at grid far-field boundary
			}
			int kD=kMax[0]-1;
			for (int mD = 0; mD <(floor((i_TOT-1)/2)+1); mD++) {
				AF[mD][kD][0][t]=Clp_nf[t]+Clf_nf[t];		//Set 1st aerodynamic force component at near field lift coeff. value at grid far-field boundary
				AF[mD][kD][1][t]=Cdp_nf[t]+Cdf_nf[t];		//Set 2nd aerodynamic force component at near field drag coeff. value at grid far-field boundary
				AF[mD][kD][2][t]=Clf_nf[t];        			//Set 3rd aerodynamic force component at near field friction lift coeff. value at grid far-field boundary
				AF[mD][kD][3][t]=Cdf_nf[t];        			//Set 4th aerodynamic force component at near field friction drag coeff. value at grid far-field boundary
			}
		}
		if (debug_mode) {
			end = chrono::steady_clock::now();
			cout << "END of including near-field aerodynamic force values in the AF data structure.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
		}
		//Compute aerodynamic force components for different sizes of the integration domain
		start = chrono::steady_clock::now();
		for (int t=0; t<Nt; t++) {
			for (int kD = k_LB; kD <k_UB; kD++) {
				for (int mD = m_LB; mD <m_UB; mD++) {
					if (debug_mode) {
						cout << "kD= "<< kD << " ; mD= " << mD << "\n";
					}
					if (mD<floor((i_TOT-1)/2)+1) {
						compute_AF(nodal_values, cell_values, aux_values, aux_values_N, flags, AF_elem, t, mD, kD, numVars, compVarsN, numAF, numAF_ADD, i_TOT, &kMax[0], NS,
									flux_scheme, grad_scheme, ParForm, limit_integr, near_ID, rho_inf, V_inf, chord, rho_id, vx_id, vz_id, X_POLE, Z_POLE, U_FORM, FT_FORM, x_CoR, z_CoR,
									ampl_x, freq_x, t0_x, ampl_z, freq_z, t0_z, ampl_a, freq_a, t0_a, MOV_GRD, time, rx_1, rz_1, rx_2, rz_2, rx_3, rz_3, rx_4, rz_4,
									dx_1, dz_1, dx_2, dz_2, dx_3, dz_3, dx_4, dz_4, a_1, a_2, a_3, a_4, WLS, AF[mD][kD]);
					}
				}
			}
		}
		if (debug_mode) {
			end = chrono::steady_clock::now();
			cout << "END of computing aerodynamic force components for different sizes of the integration domain.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
		}
	}

	//Clear data structures not used for the output phase (to leave space for tecplot output data structures to be created...)
	start = chrono::steady_clock::now();
	aux_values.clear();
	aux_values_N.clear();
	WLS.clear();
	near_ID.clear();
	if (debug_mode) {
		end = chrono::steady_clock::now();
		cout << "END of clearing data structures not used for the output phase.   (" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms)\n";
	}

	//Write out input case and settings data
	char deg=248;
	printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
	printf("%s\n",					"|                                                                                         |");
	printf("%s\n",					"|                                 UNSTEADY BREAK-FORCE 2D                                 |");
	printf("%s\n",					"|            Aerodynamic force analysis and decomposition by far-field methods            |");
	printf("%s\n",					"|                   in two-dimensional structured unsteady flow solutions                 |");
	printf("%s\n",					"|                                          v3.0                                           |");
	printf("%s\n",					"|                      Centro Italiano Ricerche Aerospaziali S.C.p.A.                     |");
	printf("%s\n",					"|                       Università degli Studi di Napoli Federico II                      |");
	printf("%s\n",					"|                                                                                         |");
	printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
	printf("%s\n",					"|                          INPUT CASE DATA and POST-PROCESSOR SETTINGS                    |");
	printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
	printf("%s%7.0d%s\n",			"|              X-VELOCITY VARIABLE ID:                                    ", vx_id, "         |");
	printf("%s%7.0d%s\n",			"|              Z-VELOCITY VARIABLE ID:                                    ", vz_id, "         |");
	if (INC) {
	printf("%s%7.0d%s\n",			"|         STATIC PRESSURE VARIABLE ID:                                    ", p_id, "         |");
	} else {
	printf("%s%7.0d%s\n",			"|            MASS DENSITY VARIABLE ID:                                    ", rho_id, "         |");
	}
	printf("%s%7.0d%s\n",			"|      STATIC TEMPERATURE VARIABLE ID:                                    ", T_id, "         |");
	printf("%s%7.0d%s\n",			"|          EDDY VISCOSITY VARIABLE ID:                                    ", mut_id, "         |");
	printf("%s%7.0d%s\n",			"|TURBULENT KINETIC ENERGY VARIABLE ID:                                    ", k_id, "         |");
	printf("%s%7.0d%s\n",			"|         GRID X-VELOCITY VARIABLE ID:                                    ", vgx_id, "         |");
	printf("%s%7.0d%s\n",			"|         GRID Z-VELOCITY VARIABLE ID:                                    ", vgz_id, "         |");
	printf("%s%7.0d%s\n",			"|  EXTERNAL ZONE-FLAGGING VARIABLE ID:                                    ", flag_id, "         |");
	printf("%s\n",					"|                                                                                         |");
	printf("%s%8.4f%s\n",			"|FREESTREAM MASS DENSITY:                                                ", rho_inf, " kg/s    |");
	printf("%s%8.4f%s%c%s\n",		"|FREESTREAM STATIC TEMPERATURE:                                          ", T_inf, " ",deg,"K      |");
	printf("%s%8.4f%s\n",			"|FREESTREAM VELOCITY MAGNITUDE:                                          ", V_inf, " m/s     |");
	printf("%s%8.4f%s\n",			"|CHORD:                                                                  ", chord, " m       |");
	printf("%s%7.2f%s\n",			"|ANGLE of ATTACK:                                                         ", AoA, " deg     |");
	printf("%s%d%s\n",			    "|EULER/NAVIER-STOKES (0/1):                                                     ", NS, "         |");
	printf("%s%d%s\n",			    "|COMPRESSIBLE/INCOMPRESSIBLE (0/1):                                             ", INC, "         |");
	printf("%s%d%s\n",			    "|INCOMPRESSIBLE FLUID or REGIME (0/1):                                          ", RHO_VAR, "         |");
	printf("%s%d%s\n",			    "|MOVING GRID (0/1):                                                             ", MOV_GRD, "         |");
	printf("%s%d%s\n",			    "|UNSTEADY FORMULATION:                                                          ", U_FORM, "         |");
	printf("%s%d%s\n",			    "|FT FORMULATION:                                                                ", FT_FORM, "         |");
	printf("%s\n",					"|                                                                                         |");
	printf("%s%7.3f%s\n",			"|                                                         F_SHOCK_TSH:    ", Fshock_tsh, "         |");
	printf("%s%7.3f%s\n",			"|                                                        F_SHOCK_TAIL:    ", Fshock_tail, "         |");
	printf("%s%7.3f%s\n",			"|                                                                K_BL:    ", K_bl, "         |");
	printf("%s%7.0f%s\n",			"|                                                           OMEGA_TSH:    ", omega_tsh, " 1/s     |");
	printf("%s%8.3f%s%c%s\n",		"|                                                            ENTR_TSH:   ", entr_tsh, " J/(kg",deg,"K)|");
	printf("%s%7.0d%s\n",			"|                                                        SHOCK_MARGIN:    ", shock_margin, "         |");
	printf("%s%7.0d%s\n",			"|                                               SHOCK_MARGIN_DRAG_TOP:    ", shock_margin_drag_kmax, "         |");
	printf("%s%7.0d%s\n",			"|                                               SHOCK_MARGIN_DRAG_MIN:    ", shock_margin_drag_mmin, "         |");
	printf("%s%7.0d%s\n",			"|                                               SHOCK_MARGIN_DRAG_MAX:    ", shock_margin_drag_mmax, "         |");
	printf("%s%8.3f%s\n",			"|                                                    SHOCK_XDS_MARGIN:   ", shock_xds_margin, " m       |");
	printf("%s%7.0d%s\n",			"|                                                         WALL_MARGIN:    ", wall_margin, "         |");
	printf("%s%7.0d%s\n",			"|                                                        BL_MARGIN_NF:    ", bl_margin_NF, "         |");
	printf("%s%7.0d%s\n",			"|                                                     BL_MARGIN_FF_UP:    ", bl_margin_FF_UP, "         |");
	printf("%s%7.0d%s\n",			"|                                                    BL_MARGIN_FF_LOW:    ", bl_margin_FF_LOW, "         |");
	printf("%s%d%s\n",				"|                                                           SMART_WDS:          ", smart_wds, "         |");
	printf("%s%d%s\n",				"|                                                         GRAD_SCHEME:          ", grad_scheme, "         |");
	printf("%s%d%s\n",				"|                                                         FLUX_SCHEME:          ", flux_scheme, "         |");
	printf("%s%d%s\n",				"|                                                           MRHO_FORM:          ", mrho_form, "         |");
	printf("%s%d%s\n",				"|                                                            PAR_FORM:          ", ParForm, "         |");
	printf("%s%d%s\n",				"|                                                           LAMB_FORM:          ", Lamb_form, "         |");
	printf("%s%d%s\n",				"|                                                          DOM_DECOMP:          ", dom_decomp, "         |");
	printf("%s%d%s\n",				"|                                                        LIMIT_INTEGR:          ", limit_integr, "         |");
	printf("%s%d%s\n",				"|                                                         THD_METHODS:          ", thd_methods, "         |");
	printf("%s%d%s\n",				"|                                                           POLE_MODE:          ", pole_mode, "         |");
	printf("%s%7.3f%s\n",			"|                                                    X_POLE (ROTATED):    ", X_POLE_IN, " m       |");
	printf("%s%7.3f%s\n",			"|                                                    Z_POLE (ROTATED):    ", Z_POLE_IN, " m       |");
	printf("%s\n",					"|                                                                                         |");
	printf("%s%7.3f%s\n",			"|                                                               X_CoR:    ", x_CoR, " m       |");
	printf("%s%7.3f%s\n",			"|                                                               Z_CoR:    ", z_CoR, " m       |");
	printf("%s%7.3f%s\n",			"|                                                   HEAVING AMPLITUDE:    ", ampl_x, " m       |");
	printf("%s%7.3f%s\n",			"|                                                   HEAVING FREQUENCY:    ", freq_x, " 1/s     |");
	printf("%s%7.3f%s\n",			"|                                                  HEAVING TIME-SHIFT:    ", t0_x, " t       |");
	printf("%s%7.3f%s\n",			"|                                                  PLUNGING AMPLITUDE:    ", ampl_z, " m       |");
	printf("%s%7.3f%s\n",			"|                                                  PLUNGING FREQUENCY:    ", freq_z, " 1/s     |");
	printf("%s%7.3f%s\n",			"|                                                 PLUNGING TIME-SHIFT:    ", t0_z, " t       |");
	printf("%s%7.3f%s\n",			"|                                                  PITCHING AMPLITUDE:    ", ampl_a, " m       |");
	printf("%s%7.3f%s\n",			"|                                                  PITCHING FREQUENCY:    ", freq_a, " 1/s     |");
	printf("%s%7.3f%s\n",			"|                                                 PITCHING TIME-SHIFT:    ", t0_a, " t       |");
	printf("%s%4d%s\n",				"|                                 NUMBER of POST-PROCESSED TIME-STEPS:       ", Nt, "         |");
	printf("%s%4d%s\n",				"|                              (ANALYSING EVERY ... INPUT TIME-STEPS):       ", every, "         |");
	printf("%s\n",					"|                                                                                         |");
	if (rescaled) {	//Input data were limited to the subset effectively used for post-processing, and m_SEL rescaled.
	printf("%s%7.0d%s\n",			"|                                                               I_SEL:    ", int(m_SEL), " rescaled|");
	} else {
	printf("%s%7.0d%s\n",			"|                                                               I_SEL:    ", int(m_SEL), "         |");
	}
	printf("%s%7.0d%s\n",			"|                                                               K_SEL:    ", int(k_SEL), "         |");
	printf("%s%7.3f%s\n",			"|                                                            (Xw_SEL):    ",
										nodal_values[int(m_SEL)][0][0][0] * cos(-AoA * M_PI / 180) + nodal_values[int(m_SEL)][0][2][0] * sin(-AoA * M_PI / 180), "         |");
	printf("%s%7.3f%s\n",			"|                                                            (Zw_SEL):    ",
										-nodal_values[int(m_SEL)][int(k_SEL)][0][0] * sin(-AoA * M_PI / 180) + nodal_values[int(m_SEL)][int(k_SEL)][2][0] * cos(-AoA * M_PI / 180), "         |");
	printf("%s%7.0d%s\n",			"|                                                             TEC_OUT:    ", TEC_OUT, "         |");
	printf("%s\n",					"|                                                                                         |");

	for (int t=0; t<Nt; t++) {
		printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
		printf("%s%4.1d%s\n",					"|                            TIME STEP #", t,"                                              |");
		//Write out NEAR-FIELD forces coefficients
		printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
		printf("%s\n",					"|                  NEAR-FIELD AERODYNAMIC FORCE (RE-COMPUTED)                             |");
		printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
		printf("%s\n",					"|                                                                                         |");
		printf("%s%7.4f%s\n",			"|TOTAL    LIFT COEFFICIENT:                                                    ", Clp_nf[t]+ Clf_nf[t], "    |");
		printf("%s%7.4f%s\n",			"|friction LIFT COEFFICIENT:                                                    ", Clf_nf[t], "    |");
		printf("%s%8.1f%s\n",			"|TOTAL    DRAG COUNTS     :                                                   ", (Cdp_nf[t] + Cdf_nf[t])*10000, "    |");
		printf("%s%8.1f%s\n",			"|friction DRAG COUNTS     :                                                   ", Cdf_nf[t]*10000, "    |");
		printf("%s\n",					"|                                                                                         |");

		//Write out FAR-FIELD forces coefficients (at specified domain size)
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		if ((k_SEL>0)&&(m_SEL>0)) {
		printf("%s\n",					"|                  FAR-FIELD AERODYNAMIC FORCE BREAK-DOWN  (VORTICITY-BASED METHOD)       |");
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		printf("%s\n",                  "|                                                                                         |");
		printf("%s%7.4f%s\n",           "|LIFT coefficient  (VORTEX FORCE):                                             ",AERO[1][t],"    |");
		printf("%s%7.4f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) :  (",AERO[numAF+9][t],")   |");
		printf("%s%7.4f%s\n",           "|                  (COMPRESSIBILITY CORRECTION - LIU FORMULATION):             ",AERO[21][t]+AERO[numAF+15][t],"    |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 1st term :  (",AERO[21][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 2nd term :  (",AERO[numAF+15][t],")   |");
		printf("%s%7.4f%s\n",           "|                  (COMPRESSIBILITY CORRECTION - ORIGINAL FORMULATION):        ",AERO[3][t],"    |");
		if (Lamb_form<2) {
		printf("%s%7.4f%s\n",           "|                  (COMPRESSIBILITY CORRECTION - MODIFIED FORMULATION):        ",AERO[5][t]+AERO[7][t]+
																														AERO[9][t]+AERO[11][t]+AERO[13][t],
																														"    |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 1st term :  (",AERO[5][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 2nd term :  (",AERO[7][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 3rd term :  (",AERO[9][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 4th term :  (",AERO[11][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 5th term :  (",AERO[13][t],")   |");
		printf("%s%7.4f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) :  (",AERO[numAF+11][t],")   |");
		printf("%s%7.4f%s\n",           "|                  (REVERSIBLE LIFT COEFFICIENT - C FORMULATION):              ",AERO[numAF+17][t]+AERO[numAF+19][t]+
																														AERO[numAF+21][t]+AERO[38][t]+
																														AERO[numAF+35][t],"    |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 1st term :  (",AERO[numAF+17][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 2nd term :  (",AERO[numAF+19][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 3rd term :  (",AERO[numAF+21][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ----------------------------------------- UNSTEADY term :  (",AERO[38][t],")   |");
		printf("%s%7.4f%s\n",           "|                  --------------------------------------------- BODY term :  (",AERO[numAF+35][t],")   |");
		}
		printf("%s%7.4f%s\n",           "|                  (LAMB VECTOR CIRC. MOMENT -  SURF. INTEGR. FORMULATION) :   ",AERO[numAF+1][t],"    |");
		printf("%s%7.4f%s\n",           "|                  (LAMB VECTOR CIRC. MOMENT - VOLUME INTEGR. FORMULATION) :   ",AERO[16][t]+AERO[19][t],"    |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 1st term :  (",AERO[16][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 2nd term :  (",AERO[19][t],")   |");
		printf("%s%7.4f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) :  (",AERO[numAF+13][t],")   |");
		if (Lamb_form<2) {
		printf("%s%7.4f%s\n",           "|                  (EXPLICIT VISCOUS STRESS TENSOR CONTRIBUTION) :             ",AERO[numAF+3][t]+AERO[numAF+5][t],"    |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 1st term :  (",AERO[numAF+3][t],")   |");
		printf("%s%7.4f%s\n",           "|                  ---------------------------------------------- 2nd term :  (",AERO[numAF+5][t],")   |");
		printf("%s%7.4f%s\n",           "|                  (ONERA COMPRESS. CORRECTION TERM - SURF. INTEGR. FORM.) :   ",-AERO[numAF+7][t],"    |");
		printf("%s%7.4f%s\n",           "|                  (ONERA COMPRESS. CORRECTION TERM -  VOL. INTEGR. FORM.) :   ",-AERO[25][t],"    |");
		printf("%s%7.4f%s\n",           "|                  (UNSTEADY TERM):                                            ",AERO[36][t],"    |");
		if (!FT_FORM) {
		printf("%s%7.4f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) :  (",AERO[numAF+33][t],")   |");
		} else {
		printf("%s%7.4f%s\n",           "|                        + Sfar surface contrib.  (inv. or mov. grid only) :  (",AERO[numAF+33][t],")   |");
		}
		printf("%s%7.4f%s\n",           "|                  (NON-INERTIAL TERM:     APPARENT FORCE from REL. MOTION):   ",AERO[40][t],"    |");
		printf("%s%7.4f%s\n",           "|                                                  + Sfar surface contrib. :   ",AERO[numAF+37][t],"    |");
		}
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		if (Lamb_form<2) {
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (LIU  FORMULATION) :                                  ",AERO[1][t]+AERO[21][t]+AERO[numAF+15][t]+
																														AERO[numAF+1][t]+
																														AERO[numAF+3][t]+AERO[numAF+5][t]+
																														AERO[numAF+9][t]+
																														AERO[36][t]+AERO[numAF+33][t]+
																														AERO[40][t]+AERO[numAF+37][t],"    |");
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (  C  FORMULATION) :                                  ",AERO[numAF+17][t]+AERO[numAF+19][t]+
																														AERO[numAF+21][t]+AERO[38][t]+
																														AERO[numAF+35][t]+AERO[numAF+1][t]+
																														AERO[numAF+3][t]+AERO[numAF+5][t]+
																														AERO[36][t]+AERO[numAF+33][t]+
																														AERO[40][t]+AERO[numAF+37][t],"    |");
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (MIX. FORMULATION: ORG. F_RHO;  VOL. INTEGR. F_S) :   ",AERO[1][t]+AERO[3][t]+
																														AERO[16][t]+AERO[19][t]+
																														AERO[numAF+3][t]+AERO[numAF+5][t]+
																														AERO[numAF+9][t]+AERO[numAF+11][t]+
																														AERO[numAF+13][t]+
																														AERO[36][t]+AERO[numAF+33][t]+
																														AERO[40][t]+AERO[numAF+37][t],"    |");
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (2014 FORMULATION: ORG. F_RHO; SURF. INTEGR. F_S) :   ",AERO[1][t]+AERO[3][t]+
																														AERO[numAF+1][t]+
																														AERO[numAF+3][t]+AERO[numAF+5][t]+
																														AERO[numAF+9][t]+AERO[numAF+11][t]+
																														AERO[36][t]+AERO[numAF+33][t]+
																														AERO[40][t]+AERO[numAF+37][t],"    |");
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (2018 FORMULATION: MOD. F_RHO;  VOL. INTEGR. F_S) :   ",AERO[1][t]+AERO[5][t]+AERO[7][t]+
																														AERO[9][t]+AERO[11][t]+AERO[13][t]+
																														AERO[16][t]+AERO[19][t]+
																														AERO[numAF+3][t]+AERO[numAF+5][t]+
																														AERO[numAF+9][t]+AERO[numAF+11][t]+
																														AERO[numAF+13][t]+
																														AERO[36][t]+AERO[numAF+33][t]+
																														AERO[40][t]+AERO[numAF+37][t],"    |");
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (2017 FORMULATION: MOD. F_RHO; SURF. INTEGR. F_S) :   ",AERO[1][t]+AERO[5][t]+AERO[7][t]+
																														AERO[9][t]+AERO[11][t]+AERO[13][t]+
																														AERO[numAF+1][t]+AERO[numAF+3][t]+AERO[numAF+5][t]+
																														AERO[numAF+9][t]+AERO[numAF+11][t]+
																														AERO[36][t]+AERO[numAF+33][t]+
																														AERO[40][t]+AERO[numAF+37][t],"    |");
/* 		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (2021 FORMULATION: ONERA     ; SURF. INTEGR. F_S) :   ",AERO[1][t]+AERO[5][t]+AERO[7][t]+
																														AERO[9][t]+AERO[11][t]+AERO[13][t]-
																														AERO[numAF+7][t]+
																														AERO[numAF+9][t]+AERO[numAF+11][t]+
																														AERO[36][t]+AERO[numAF+33][t],"    |"); */
		} else {
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (LIU  FORMULATION) :                                  ",AERO[1][t]+AERO[21][t]+AERO[numAF+15][t]+
																														AERO[numAF+1][t]+
																														AERO[numAF+9][t]+Clf_nf[t],"    |");
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (MIX. FORMULATION: ORG. F_RHO;  VOL. INTEGR. F_S) :   ",AERO[1][t]+AERO[3][t]+
																														AERO[16][t]+AERO[19][t]+
																														AERO[numAF+9][t]+AERO[numAF+11][t]+
																														AERO[numAF+13][t]+Clf_nf[t],"    |");
		printf("%s%7.4f%s\n",           "|TOTAL LIFT COEFFICIENT  (2014 FORMULATION: ORG. F_RHO; SURF. INTEGR. F_S) :   ",AERO[1][t]+AERO[3][t]+
																														AERO[numAF+1][t]+																														
																														AERO[numAF+9][t]+AERO[numAF+11][t]+
																														Clf_nf[t],"    |");
		}
 		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		printf("%s\n",                  "|                                                                                         |");
		printf("%s%8.1f%s\n",           "|DRAG counts       (VORTEX FORCE):                                            ",10000*(AERO[0][t]),"    |");
		printf("%s%8.1f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) : (",10000*(AERO[numAF+8][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  (COMPRESSIBILITY CORRECTION - LIU FORMULATION):            ",10000*(AERO[20][t]+AERO[numAF+14][t]),"    |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 1st term : (",10000*(AERO[20][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 2nd term : (",10000*(AERO[numAF+14][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  (COMPRESSIBILITY CORRECTION - ORIGINAL FORMULATION):       ",10000*(AERO[2][t]),"    |");
		if (Lamb_form<2) {
		printf("%s%8.1f%s\n",           "|                  (COMPRESSIBILITY CORRECTION - MODIFIED FORMULATION):       ",10000*(AERO[4][t]+AERO[6][t]+
																														AERO[8][t]+AERO[10][t]+AERO[12][t]),
																														"    |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 1st term : (",10000*(AERO[4][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 2nd term : (",10000*(AERO[6][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 3rd term : (",10000*(AERO[8][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 4th term : (",10000*(AERO[10][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 5th term : (",10000*(AERO[12][t]),")   |");
		printf("%s%8.1f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) : (",10000*(AERO[numAF+10][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  (REVERSIBLE DRAG COUNTS - C FORMULATION):                  ",10000*(AERO[numAF+16][t]+AERO[numAF+18][t]+
																														AERO[numAF+20][t]+AERO[37][t]+
																														AERO[numAF+34][t]),"    |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 1st term : (",10000*(AERO[numAF+16][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 2nd term : (",10000*(AERO[numAF+18][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 3rd term : (",10000*(AERO[numAF+20][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------- (All terms, VISCOUS region only) : (",10000*(AERO[numAF+26][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ------------------------- (All terms, WAVE region only) : (",10000*(AERO[numAF+42][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ----------------------------------------- UNSTEADY term : (",10000*(AERO[37][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  --------------------------------------------- BODY term : (",10000*(AERO[numAF+34][t]),")   |");
		}
		printf("%s%8.1f%s\n",           "|                  (LAMB VECTOR CIRC. MOMENT - SURF. INTEGR.  FORMULATION) :  ",10000*(AERO[numAF][t]),"    |");
		printf("%s%8.1f%s\n",			"|                  --------------------------------- (VISCOUS region only) :  ",10000*(AERO[numAF+22][t]),"    |");
		printf("%s%8.1f%s\n",			"|                  --------------------------------- (VISCOUS  WAKE  only) :  ",10000*(AERO[numAF+28][t]),"    |");
		printf("%s%8.1f%s\n",			"|                  ------------------------------------ (WAVE region only) :  ",10000*(AERO[numAF+38][t]),"    |");
		printf("%s%8.1f%s\n",           "|                  (LAMB VECTOR CIRC. MOMENT - VOLUME INTEGR. FORMULATION) :  ",10000*(AERO[14][t]+AERO[15][t]+AERO[33][t]+
																														AERO[17][t]+AERO[18][t]+AERO[34][t]),"    |");
		printf("%s%8.1f%s\n",           "|                  --------------------------------- 1st term ;       WAVE : (",10000*(AERO[14][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  --------------------------------- 1st term ;    VISCOUS : (",10000*(AERO[15][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  --------------------------------- 1st term ;   SPURIOUS : (",10000*(AERO[33][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  --------------------------------- 2nd term ;       WAVE : (",10000*(AERO[17][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  --------------------------------- 2nd term ;    VISCOUS : (",10000*(AERO[18][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  --------------------------------- 2nd term ;   SPURIOUS : (",10000*(AERO[34][t]),")   |");
		printf("%s%8.1f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) : (",10000*(AERO[numAF+12][t]),")   |");
		if (Lamb_form<2) {
		printf("%s%8.1f%s\n",           "|                  (EXPLICIT VISCOUS STRESS TENSOR CONTRIBUTION):             ",10000*(AERO[numAF+2][t]+AERO[numAF+4][t]),"    |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 1st term : (",10000*(AERO[numAF+2][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  ---------------------------------------------- 2nd term : (",10000*(AERO[numAF+4][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  (ONERA COMPRESS. CORR. TERM - SURF. INTEGR. FORM.)   :  +/-",10000*(AERO[numAF+6][t]),"    |");
		printf("%s%8.1f%s\n",			"|                  --------------------------------- (VISCOUS region only) :  ",10000*(AERO[numAF+24][t]),"    |");
		printf("%s%8.1f%s\n",			"|                  --------------------------------- (VISCOUS  WAKE  only) :  ",10000*(AERO[numAF+30][t]),"    |");
		printf("%s%8.1f%s\n",			"|                  ------------------------------------ (WAVE region only) :  ",10000*(AERO[numAF+40][t]),"    |");
		//printf("%s%8.1f%s\n",			"|                  ------------------------------------ (WAVE region only) :  ",10000*(AERO[32]),"    |");
		printf("%s%8.1f%s\n",           "|                  (ONERA COMPRESS. CORR. TERM -  VOL. INTEGR. FORM.)   :  +/-",10000*(AERO[22][t]+
																																AERO[23][t]+AERO[24][t]),"    |");
		printf("%s%8.1f%s\n",           "|                  -------------------------------- all terms ;       WAVE : (",10000*(AERO[22][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  -------------------------------- all terms ;    VISCOUS : (",10000*(AERO[23][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  -------------------------------- all terms ;   SPURIOUS : (",10000*(AERO[24][t]),")   |");
		printf("%s%8.1f%s\n",           "|                  (UNSTEADY TERM):                                           ",10000*(AERO[35][t]),"    |");
		if (!FT_FORM) {
		printf("%s%8.1f%s\n",           "|                        + Body surface contrib.  (inv. or mov. grid only) : (",10000*(AERO[numAF+32][t]),")   |");
		} else {
		printf("%s%8.1f%s\n",           "|                        + Sfar surface contrib.  (inv. or mov. grid only) : (",10000*(AERO[numAF+32][t]),")   |");
		}
		printf("%s%8.1f%s\n",           "|                  (NON-INERTIAL TERM:     APPARENT FORCE from REL. MOTION):  ",10000*(AERO[39][t]),"    |");
		printf("%s%8.1f%s\n",           "|                                                  + Sfar surface contrib. :  ",10000*(AERO[numAF+36][t]),"    |");
		}
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		if (Nt==1) {
		if (Lamb_form<2) {
		if (ONERA_corr) {
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (     LIU FORMULATION)  ------------------------ :  ",10000*(AERO[0][t]+AERO[20][t]+AERO[numAF+14][t]+
																														AERO[numAF+8][t]-AERO[numAF+6][t]),"    |");
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (       C FORMULATION)  ------------------------ :  ",10000*(AERO[numAF+16][t]+AERO[numAF+18][t]+
																														AERO[numAF+20][t]-AERO[numAF+6][t]),"    |");
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (ORIGINAL FORMULATION)  ------------------------ :  ",10000*(AERO[0][t]+AERO[2][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														-AERO[numAF+6][t]),"    |");
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (MODIFIED FORMULATION)  ------------------------ :  ",10000*(AERO[0][t]+AERO[4][t]+AERO[6][t]+
																														AERO[8][t]+AERO[10][t]+AERO[12][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														-AERO[numAF+6][t]),"    |");
		} else {
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (     LIU FORMULATION)  ------------------------ :  ",10000*(AERO[0][t]+AERO[20][t]+AERO[numAF+14][t]+
																														AERO[numAF+8][t]),"    |");
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (       C FORMULATION)  ------------------------ :  ",10000*(AERO[numAF+16][t]+AERO[numAF+18][t]+
																														AERO[numAF+20][t]),"    |");
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (ORIGINAL FORMULATION)  ------------------------ :  ",10000*(AERO[0][t]+AERO[2][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]),"    |");
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (MODIFIED FORMULATION)  ------------------------ :  ",10000*(AERO[0][t]+AERO[4][t]+AERO[6][t]+
																														AERO[8][t]+AERO[10][t]+AERO[12][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]),"    |");
		printf("%s%8.1f%s\n",           "|INDUCED DRAG COUNTS      (ONERA    FORMULATION)  ------------------------ :  ",10000*(AERO[0][t]+AERO[4][t]+AERO[6][t]+
																														AERO[8][t]+AERO[10][t]+AERO[12][t]-
																														AERO[numAF+6][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]),"    |");
		}
		} else {
			//To be implemented
		}
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		if (Lamb_form<2) {
		if (ONERA_corr) {
		printf("%s%8.1f%s\n",           "|PARASS. DRAG COUNTS      (SURF. INTEGR.  FORMULATION) ------------- TOTAL :  ",10000*(AERO[numAF][t]+AERO[numAF+6][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]),"    |");
		printf("%s%8.1f%s\n",			"|                         ------------------------------------------- WAVE :  ",10000*(AERO[numAF+38][t]+AERO[numAF+40][t]),"    |");
		printf("%s%8.1f%s\n",			"|                         ---------------------------------------- VISCOUS :  ",10000*(AERO[numAF+22][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]+
																																AERO[numAF+24][t]),"    |");
		printf("%s%8.1f%s\n",			"|                         --------------------------------------- SPURIOUS :  ",10000*(AERO[numAF][t]+AERO[numAF+6][t]-AERO[numAF+22][t]-
																																AERO[numAF+24][t]-AERO[numAF+38][t]-AERO[numAF+40][t]),"    |");
		printf("%s%8.1f%s\n",           "|PARASS. DRAG COUNTS      (VOLUME INTEGR. FORMULATION) ------------- TOTAL :  ",10000*(AERO[14][t]+AERO[15][t]+AERO[33][t]+
																																AERO[17][t]+AERO[18][t]+AERO[34][t]+
																																AERO[22][t]+AERO[23][t]+AERO[24][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]+
																																AERO[numAF+12][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         ------------------------------------------- WAVE :  ",10000*(AERO[14][t]+AERO[17][t]+
																																AERO[22][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         ---------------------------------------- VISCOUS :  ",10000*(AERO[15][t]+AERO[18][t]+AERO[23][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         --------------------------------------- SPURIOUS :  ",10000*(AERO[33][t]+AERO[34][t]+AERO[24][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         ----------------- BODY  (inv. or mov. grid only) :  ",10000*(AERO[numAF+12][t]),"    |");
		} else {
		printf("%s%8.1f%s\n",           "|PARASS. DRAG COUNTS      (SURF. INTEGR.  FORMULATION) ------------- TOTAL :  ",10000*(AERO[numAF][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]),"    |");
		printf("%s%8.1f%s\n",			"|                         ------------------------------------------- WAVE :  ",10000*(AERO[numAF+38][t]),"    |");
		printf("%s%8.1f%s\n",			"|                         ---------------------------------------- VISCOUS :  ",10000*(AERO[numAF+22][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]),"    |");
		printf("%s%8.1f%s\n",			"|                         --------------------------------------- SPURIOUS :  ",10000*(AERO[numAF][t]-AERO[numAF+22][t]-
																																-AERO[numAF+38][t]),"    |");
		printf("%s%8.1f%s\n",           "|PARASS. DRAG COUNTS      (VOLUME INTEGR. FORMULATION) ------------- TOTAL :  ",10000*(AERO[14][t]+AERO[15][t]+AERO[33][t]+
																																AERO[17][t]+AERO[18][t]+AERO[34][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]+
																																AERO[numAF+12][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         ------------------------------------------- WAVE :  ",10000*(AERO[14][t]+AERO[17][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         ---------------------------------------- VISCOUS :  ",10000*(AERO[15][t]+AERO[18][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         --------------------------------------- SPURIOUS :  ",10000*(AERO[33][t]+AERO[34][t]),"    |");
		printf("%s%8.1f%s\n",           "|                         ----------------- BODY  (inv. or mov. grid only) :  ",10000*(AERO[numAF+12][t]),"    |");
		printf("%s%8.1f%s\n",           "|PARASS. DRAG COUNTS      (ONERA SURF. INTEGR. FORMULATION ) ------- TOTAL :  ",10000*(AERO[numAF][t]+
																																AERO[numAF+2][t]+AERO[numAF+4][t]+
																																AERO[numAF+6][t]),"    |");
		}
		} else {
			//To be implemented
		}
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		}
		if (Lamb_form<2) {
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (LIU  FORMULATION):                                 ",10000*(AERO[0][t]+AERO[20][t]+
																														AERO[numAF+14][t]+AERO[numAF][t]+
																														AERO[numAF+2][t]+AERO[numAF+4][t]+
																														AERO[numAF+8][t]+
																														AERO[35][t]+AERO[numAF+32][t]+
																														AERO[39][t]+AERO[numAF+36][t]),"    |");
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (  C  FORMULATION):                                 ",10000*(AERO[numAF+16][t]+AERO[numAF+18][t]+
																														AERO[numAF+20][t]+AERO[37][t]+
																														AERO[numAF+34][t]+AERO[numAF][t]+
																														AERO[numAF+2][t]+AERO[numAF+4][t]+
																														AERO[35][t]+AERO[numAF+32][t]+
																														AERO[39][t]+AERO[numAF+36][t]),"    |");
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (MIX. FORMULATION: ORG. F_RHO;  VOL. INTEGR. F_S):  ",10000*(AERO[0][t]+AERO[2][t]+AERO[14][t]+
																														AERO[15][t]+AERO[17][t]+AERO[18][t]+
																														AERO[33][t]+AERO[34][t]+
																														AERO[numAF+2][t]+AERO[numAF+4][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														AERO[numAF+12][t]+
																														AERO[35][t]+AERO[numAF+32][t]+
																														AERO[39][t]+AERO[numAF+36][t]),"    |");
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (2014 FORMULATION: ORG. F_RHO; SURF. INTEGR. F_S):  ",10000*(AERO[0][t]+AERO[2][t]+
																														AERO[numAF][t]+
																														AERO[numAF+2][t]+AERO[numAF+4][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														AERO[35][t]+AERO[numAF+32][t]+
																														AERO[39][t]+AERO[numAF+36][t]),"    |");
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (2018 FORMULATION: MOD. F_RHO;  VOL. INTEGR. F_S):  ",10000*(AERO[0][t]+AERO[4][t]+AERO[6][t]+
																														AERO[8][t]+AERO[10][t]+AERO[12][t]+
																														AERO[14][t]+AERO[15][t]+AERO[17][t]+
																														AERO[18][t]+AERO[33][t]+AERO[34][t]+
																														AERO[numAF+2][t]+AERO[numAF+4][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														AERO[numAF+12][t]+
																														AERO[35][t]+AERO[numAF+32][t]+
																														AERO[39][t]+AERO[numAF+36][t]),"    |");
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (2017 FORMULATION: MOD. F_RHO; SURF. INTEGR. F_S):  ",10000*(AERO[0][t]+AERO[4][t]+AERO[6][t]+
																														AERO[8][t]+AERO[10][t]+AERO[12][t]+
																														AERO[numAF][t]+
																														AERO[numAF+2][t]+AERO[numAF+4][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														AERO[35][t]+AERO[numAF+32][t]+
																														AERO[39][t]+AERO[numAF+36][t]),"    |");
		} else {
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (LIU  FORMULATION):                                 ",10000*(AERO[0][t]+AERO[20][t]+
																														AERO[numAF+14][t]+AERO[numAF][t]+																														
																														AERO[numAF+8][t]+Cdf_nf[t]),"    |");
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (MIX. FORMULATION: ORG. F_RHO;  VOL. INTEGR. F_S):  ",10000*(AERO[0][t]+AERO[2][t]+AERO[14][t]+
																														AERO[15][t]+AERO[17][t]+AERO[18][t]+
																														AERO[33][t]+AERO[34][t]+																														
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														AERO[numAF+12][t]+Cdf_nf[t]),"    |");
		printf("%s%8.1f%s\n",           "|  TOTAL DRAG COUNTS      (2014 FORMULATION: ORG. F_RHO; SURF. INTEGR. F_S):  ",10000*(AERO[0][t]+AERO[2][t]+
																														AERO[numAF][t]+
																														AERO[numAF+8][t]+AERO[numAF+10][t]+
																														Cdf_nf[t]),"    |");
		}
		printf("%s\n",					"|                                                                                         |");
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		if (Nt==1) {
		printf("%s\n",					"|                  FAR-FIELD AERODYNAMIC DRAG BREAK-DOWN  (THERMODYNAMIC METHOD)          |");
		printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
		printf("%s\n",					"|                                                                                         |");
		printf("%s%8.1f%s\n",			"|  TOTAL DRAG COUNTS       (DESTARAC - VAN DER VOOREN  formulation)  ----- :  ", 10000 * (AERO[26][t] + AERO[27][t] + 
																																AERO[28][t]), "    |");
		printf("%s%8.1f%s\n",			"|                          ------------------------------------------ WAVE :  ", 10000 * (AERO[26][t]), "    |");
		printf("%s%8.1f%s\n",			"|                          --------------------------------------- VISCOUS :  ", 10000 * (AERO[27][t]), "    |");
		printf("%s%8.1f%s\n",			"|                          -------------------------------------- SPURIOUS :  ", 10000 * (AERO[28][t]), "    |");
		printf("%s%8.1f%s\n",			"|  TOTAL DRAG COUNTS       (PAPARONE - TOGNACCINI      formulation)  ------:  ", 10000 * (AERO[29][t] + AERO[30][t] +
																																AERO[31][t]), "    |");
		printf("%s%8.1f%s\n",			"|                          ------------------------------------------ WAVE :  ", 10000 * (AERO[29][t]), "    |");
		printf("%s%8.1f%s\n",			"|                          --------------------------------------- VISCOUS :  ", 10000 * (AERO[30][t]), "    |");
		printf("%s%8.1f%s\n",			"|                          -------------------------------------- SPURIOUS :  ", 10000 * (AERO[31][t]), "    |");
		printf("%s\n",					"|                                                                                         |");
		printf("%s\n",					"|-----------------------------------------------------------------------------------------|");
		}
		} else {
		printf("%s\n",                  "| FAR-FIELD AERODYNAMIC FORCE BREAK-DOWN is being written in parametric Tecplot format... |");
		printf("%s\n",                  "|-----------------------------------------------------------------------------------------|");
		}
	}

	//Time history
	FILE* histFile;
	if ((Nt>=3)&&(k_SEL>0)&&(m_SEL>0)) {	//Only written in case of unsteady post-processing with a single integration domain.
											//Tecplot output will contain similar data in case of variable integration domain.
		string HistFileName=argv[2];
                HistFileName.resize(HistFileName.length()-3);
                HistFileName=HistFileName+"dat";
		histFile = fopen (HistFileName.c_str(),"w");
		if (!FT_FORM) {
			fprintf (histFile, "%s%s%s%s%s\n",
				"t/t_adv.\tCl_nf\tCd_nf\tCl_nf_fric\tCd_nf_fric\tCl_ff_LIU\tCd_ff_LIU\tCl_ff_2014\tCd_ff_2014\tCl_ff_MIX\tCd_ff_MIX\tCl_ff_C\tCd_ff_C\tCl_ff_2017\tCd_ff_2017\tCl_ff_2018\tCd_ff_2018\t",
				"Cl_LAMB\tCl_RHO\tCl_M-RHO\tCl_M-RHO-C\tCl_RHO+LAMB\tCl_S\tCl_S-OMEGA\tCl_MU-T\tCl_GRAD-RHO\tCl_GRAD-RHO-OMEGA\tCl_UNS\tCl_APP\t",
				"Cd_LAMB\tCd_RHO\tCd_M-RHO\tCd_M-RHO-C\tCd_RHO+LAMB\tCd_S\tCd_S-OMEGA\tCd_MU-T\tCd_GRAD-RHO\tCd_GRAD-RHO-OMEGA\tCd_UNS\tCd_APP\t",
				"Cl_LAMB-BODY\tCl_COMPR-BODY\tCl_RHO+LAMB-BODY\tCl_S-BODY\tCl_UNS-BODY\tCl_APP-SFAR\t",
				"Cd_LAMB-BODY\tCd_COMPR-BODY\tCd_RHO+LAMB-BODY\tCd_S-BODY\tCd_UNS-BODY\tCd_APP-SFAR");
		} else {
			fprintf (histFile, "%s%s%s%s%s\n",
				"t/t_adv.\tCl_nf\tCd_nf\tCl_nf_fric\tCd_nf_fric\tCl_ff_LIU\tCd_ff_LIU\tCl_ff_2014\tCd_ff_2014\tCl_ff_MIX\tCd_ff_MIX\tCl_ff_C\tCd_ff_C\tCl_ff_2017\tCd_ff_2017\tCl_ff_2018\tCd_ff_2018\t",
				"Cl_LAMB\tCl_RHO\tCl_M-RHO\tCl_M-RHO-C\tCl_RHO+LAMB\tCl_S\tCl_S-OMEGA\tCl_MU-T\tCl_GRAD-RHO\tCl_GRAD-RHO-OMEGA\tCl_UNS\tCl_APP\t",
				"Cd_LAMB\tCd_RHO\tCd_M-RHO\tCd_M-RHO-C\tCd_RHO+LAMB\tCd_S\tCd_S-OMEGA\tCd_MU-T\tCd_GRAD-RHO\tCd_GRAD-RHO-OMEGA\tCd_UNS\tCd_APP\t",
				"Cl_LAMB-BODY\tCl_COMPR-BODY\tCl_RHO+LAMB-BODY\tCl_S-BODY\tCl_UNS-SFAR\tCl_APP-SFAR\t",
				"Cd_LAMB-BODY\tCd_COMPR-BODY\tCd_RHO+LAMB-BODY\tCd_S-BODY\tCd_UNS-SFAR\tCd_APP-SFAR");
		}
		for (int t=0; t<Nt; t++) {
			fprintf (histFile, "%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t",
				time[t]*V_inf/chord, Clp_nf[t]+ Clf_nf[t], Cdp_nf[t] + Cdf_nf[t], Clf_nf[t], Cdf_nf[t]);
			if (Lamb_form<2) {
				fprintf (histFile, "%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t",
					//LIU
					AERO[1][t]+AERO[21][t]+AERO[numAF+15][t]+AERO[numAF+1][t]+AERO[numAF+3][t]+AERO[numAF+5][t]+AERO[numAF+9][t]+AERO[36][t]+AERO[numAF+33][t]+AERO[40][t]+AERO[numAF+37][t],
					AERO[0][t]+AERO[20][t]+AERO[numAF+14][t]+AERO[numAF][t]+AERO[numAF+2][t]+AERO[numAF+4][t]+AERO[numAF+8][t]+AERO[35][t]+AERO[numAF+32][t]+AERO[39][t]+AERO[numAF+36][t],
					//2014
					AERO[1][t]+AERO[3][t]+AERO[numAF+1][t]+AERO[numAF+3][t]+AERO[numAF+5][t]+AERO[numAF+9][t]+AERO[numAF+11][t]+AERO[36][t]+AERO[numAF+33][t]+AERO[40][t]+AERO[numAF+37][t],
					AERO[0][t]+AERO[2][t]+AERO[numAF][t]+AERO[numAF+2][t]+AERO[numAF+4][t]+AERO[numAF+8][t]+AERO[numAF+10][t]+AERO[35][t]+AERO[numAF+32][t]+AERO[39][t]+AERO[numAF+36][t],
					//MIX
					AERO[1][t]+AERO[3][t]+AERO[16][t]+AERO[19][t]+AERO[numAF+3][t]+AERO[numAF+5][t]+AERO[numAF+9][t]+AERO[numAF+11][t]+AERO[numAF+13][t]+AERO[36][t]+AERO[numAF+33][t]+AERO[40][t]+AERO[numAF+37][t],
					AERO[0][t]+AERO[2][t]+AERO[14][t]+AERO[15][t]+AERO[17][t]+AERO[18][t]+AERO[33][t]+AERO[34][t]+AERO[numAF+2][t]+AERO[numAF+4][t]+AERO[numAF+8][t]+AERO[numAF+10][t]+AERO[numAF+12][t]+AERO[35][t]+AERO[numAF+32][t]+AERO[39][t]+AERO[numAF+36][t],
					//C
					AERO[numAF+17][t]+AERO[numAF+19][t]+AERO[numAF+21][t]+AERO[38][t]+AERO[numAF+35][t]+AERO[numAF+1][t]+AERO[numAF+3][t]+AERO[numAF+5][t]+AERO[36][t]+AERO[numAF+33][t]+AERO[40][t]+AERO[numAF+37][t],
					AERO[numAF+16][t]+AERO[numAF+18][t]+AERO[numAF+20][t]+AERO[37][t]+AERO[numAF+34][t]+AERO[numAF][t]+AERO[numAF+2][t]+AERO[numAF+4][t]+AERO[35][t]+AERO[numAF+32][t]+AERO[39][t]+AERO[numAF+36][t],
					//2017
					AERO[1][t]+AERO[5][t]+AERO[7][t]+AERO[9][t]+AERO[11][t]+AERO[13][t]+AERO[numAF+1][t]+AERO[numAF+3][t]+AERO[numAF+5][t]+AERO[numAF+9][t]+AERO[numAF+11][t]+AERO[36][t]+AERO[numAF+33][t]+AERO[40][t]+AERO[numAF+37][t],
					AERO[0][t]+AERO[4][t]+AERO[6][t]+AERO[8][t]+AERO[10][t]+AERO[12][t]+AERO[numAF][t]+AERO[numAF+2][t]+AERO[numAF+4][t]+AERO[numAF+8][t]+AERO[numAF+10][t]+AERO[35][t]+AERO[numAF+32][t]+AERO[39][t]+AERO[numAF+36][t],
					//2018
					AERO[1][t]+AERO[5][t]+AERO[7][t]+AERO[9][t]+AERO[11][t]+AERO[13][t]+AERO[16][t]+AERO[19][t]+AERO[numAF+3][t]+AERO[numAF+5][t]+AERO[numAF+9][t]+AERO[numAF+11][t]+AERO[numAF+13][t]+AERO[36][t]+AERO[numAF+33][t]+AERO[40][t]+AERO[numAF+37][t],
					AERO[0][t]+AERO[4][t]+AERO[6][t]+AERO[8][t]+AERO[10][t]+AERO[12][t]+AERO[14][t]+AERO[15][t]+AERO[17][t]+AERO[18][t]+AERO[33][t]+AERO[34][t]+AERO[numAF+2][t]+AERO[numAF+4][t]+AERO[numAF+8][t]+AERO[numAF+10][t]+AERO[numAF+12][t]+AERO[35][t]+AERO[numAF+32][t]+AERO[39][t]+AERO[numAF+36][t]);
			} else {
				fprintf (histFile, "%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t",
					//LIU
					AERO[1][t]+AERO[21][t]+AERO[numAF+15][t]+AERO[numAF+1][t]+AERO[numAF+9][t]+Clf_nf[t],
					AERO[0][t]+AERO[20][t]+AERO[numAF+14][t]+AERO[numAF][t]+AERO[numAF+8][t]+Cdf_nf[t],
					//2014
					AERO[1][t]+AERO[3][t]+AERO[numAF+1][t]+AERO[numAF+9][t]+AERO[numAF+11][t]+Clf_nf[t],
					AERO[0][t]+AERO[2][t]+AERO[numAF][t]+AERO[numAF+8][t]+AERO[numAF+10][t]+Cdf_nf[t],
					//MIX
					AERO[1][t]+AERO[3][t]+AERO[16][t]+AERO[19][t]+AERO[numAF+9][t]+AERO[numAF+11][t]+AERO[numAF+13][t]+Clf_nf[t],
					AERO[0][t]+AERO[2][t]+AERO[14][t]+AERO[15][t]+AERO[17][t]+AERO[18][t]+AERO[33][t]+AERO[34][t]+AERO[numAF+8][t]+AERO[numAF+10][t]+AERO[numAF+12][t]+Cdf_nf[t],
					//C
					NAN,
					NAN,
					//2017
					NAN,
					NAN,
					//2018
					NAN,
					NAN);
			}
			fprintf (histFile, "%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t",
				AERO[1][t], AERO[21][t]+AERO[numAF+15][t], AERO[3][t], AERO[5][t]+AERO[7][t]+AERO[9][t]+AERO[11][t]+AERO[13][t], AERO[numAF+17][t]+AERO[numAF+19][t]+AERO[numAF+21][t]+AERO[38][t],
					AERO[numAF+1][t], AERO[16][t]+AERO[19][t], AERO[numAF+3][t]+AERO[numAF+5][t], -AERO[numAF+7][t], -AERO[25][t], AERO[36][t], AERO[40][t]);
			fprintf (histFile, "%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t",
				AERO[0][t], AERO[20][t]+AERO[numAF+14][t], AERO[2][t], AERO[4][t]+AERO[6][t]+AERO[8][t]+AERO[10][t]+AERO[12][t], AERO[numAF+16][t]+AERO[numAF+18][t]+AERO[numAF+20][t]+AERO[37][t],
					AERO[numAF][t], AERO[14][t]+AERO[15][t]+AERO[33][t]+AERO[17][t]+AERO[18][t]+AERO[34][t], AERO[numAF+2][t]+AERO[numAF+4][t], AERO[numAF+6][t], AERO[22][t]+AERO[23][t]+AERO[24][t],
					AERO[35][t], AERO[39][t]);
			fprintf (histFile, "%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t",
				AERO[numAF+9][t], AERO[numAF+11][t], AERO[numAF+35][t], AERO[numAF+13][t], AERO[numAF+33][t], AERO[numAF+37][t]);
			fprintf (histFile, "%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E\t%22.15E",
				AERO[numAF+8][t], AERO[numAF+10][t], AERO[numAF+34][t], AERO[numAF+12][t], AERO[numAF+32][t], AERO[numAF+36][t]);
			if (t<Nt-1) {
				fprintf (histFile, "\n");
			}
		}
		fclose (histFile);
	}

	//TECPLOT output
	if (TEC_OUT>0) {
		string out_var_names[compVarsN + compVarsC + numAF + numFlags] = {
			"Static Pressure",
			"Specific kinetic energy per unit mass",
			"Cell area",
			"Dynamic (molecular) viscosity",
			"Velocity gradient (du/dx)",
			"Velocity gradient (du/dz)",
			"Velocity gradient (dw/dx)",
			"Velocity gradient (dw/dz)",
			"Static pressure gradient (x-component)",
			"Static pressure gradient (z-component)",
			"Mass density gradient (x-component)",
			"Mass density gradient (z-component)",
			"Specific kinetic energy gradient (x-component)",
			"Specific kinetic energy gradient (z-component)",
			"Vorticity (y-component)",
			"Lamb vector (x-component)",
			"Lamb vector (z-component)",
			"Cell-centroid (x-coordinate)",
			"Cell-centroid (z-coordinate)",
			"Shock function",
			"Eddy viscosity BL sensor",
			"Turbulence dissipation frequency (omega)",
			"Entropy variation (w.r.t. freestream)",
			"Elem. vortex force coeff. (x-component)",
			"Elem. vortex force coeff. (z-component)",
			"Elem. compr. corr. coeff. (x-component)",
			"Elem. compr. corr. coeff. (z-component)",
			"Elem. modif. compr. corr. coeff. (1st term | x-component)",
			"Elem. modif. compr. corr. coeff. (1st term | z-component)",
			"Elem. modif. compr. corr. coeff. (2nd term | x-component)",
			"Elem. modif. compr. corr. coeff. (2nd term | z-component)",
			"Elem. modif. compr. corr. coeff. (3rd term | x-component)",
			"Elem. modif. compr. corr. coeff. (3rd term | z-component)",
			"Elem. modif. compr. corr. coeff. (4th term | x-component)",
			"Elem. modif. compr. corr. coeff. (4th term | z-component)",
			"Elem. modif. compr. corr. coeff. (5th term | x-component)",
			"Elem. modif. compr. corr. coeff. (5th term | z-component)",
			"Elem. Lamb vect. circ. mom. coeff. (1st term | SW region | x-component)",
			"Elem. Lamb vect. circ. mom. coeff. (1st term | VISCOUS region | x-component)",
			"Elem. Lamb vect. circ. mom. coeff. (1st term | z-component)",
			"Elem. Lamb vect. circ. mom. coeff. (2nd term | SW region | x-component)",
			"Elem. Lamb vect. circ. mom. coeff. (2nd term | VISCOUS region | x-component)",
			"Elem. Lamb vect. circ. mom. coeff. (2nd term | z-component)",
			"Elem. Liu Compress. corr. coeff. (1st term | x-component)",
			"Elem. Liu Compress. corr. coeff. (1st term | z-component)",
			"Elem. ONERA Compress. corr. coeff. (Volume formul. | WAVE region | x-component)",
			"Elem. ONERA Compress. corr. coeff. (Volume formul. | VISCOUS region | x-component)",
			"Elem. ONERA Compress. corr. coeff. (Volume formul. | SPURIOUS region | x-component)",
			"Elem. ONERA Compress. corr. coeff. (Volume formul. | z-component)",
			"Elem. drag coeff. (TD method 1 - WAVE)",
			"Elem. drag coeff. (TD method 1 - VISCOUS)",
			"Elem. drag coeff. (TD method 1 - SPURIOUS)",
			"Elem. entropy drag coeff. (TD method 1 - WAVE)",
			"Elem. entropy drag coeff. (TD method 1 - VISCOUS)",
			"Elem. entropy drag coeff. (TD method 1 - SPURIOUS)",
			"Elem. Lamb vect. circ. mom. coeff. (WAVE region surface integral)",
			"Elem. Lamb vect. circ. mom. coeff. (1st term | SPURIOUS region | x-component)",
			"Elem. Lamb vect. circ. mom. coeff. (2nd term | SPURIOUS region | x-component)",
			"Elem. Unsteady term coeff. (volume integrand | x-component)",
			"Elem. Unsteady term coeff. (volume integrand | z-component)",
			"Elem. Unsteady term coeff. (C-formulation volume integrand | x-component)",
			"Elem. Unsteady term coeff. (C-formulation volume integrand | z-component)",
			"Elem. Non-Inertial term coeff. (volume integrand | x-component)",
			"Elem. Non-Inertial term coeff. (volume integrand | z-component)",
			"Shock wave domain FLAG",
			"Viscous domain FLAG",
			"Wave drag domain FLAG (vortical)",
			"Wave drag domain FLAG (thermod.)"
		};
		string AF_var_names[numAF + numAF_ADD] = {
			"Vortex force coeff. (x-component)",
			"Vortex force coeff. (z-component)",
			"Compr. corr. coeff. (x-component)",
			"Compr. corr. coeff. (z-component)",
			"Modif. compr. corr. coeff. (1st term | x-component)",
			"Modif. compr. corr. coeff. (1st term | z-component)",
			"Modif. compr. corr. coeff. (2nd term | x-component)",
			"Modif. compr. corr. coeff. (2nd term | z-component)",
			"Modif. compr. corr. coeff. (3rd term | x-component)",
			"Modif. compr. corr. coeff. (3rd term | z-component)",
			"Modif. compr. corr. coeff. (4th term | x-component)",
			"Modif. compr. corr. coeff. (4th term | z-component)",
			"Modif. compr. corr. coeff. (5th term | x-component)",
			"Modif. compr. corr. coeff. (5th term | z-component)",
			"Lamb vect. circ. mom. coeff. VOLUME form. (1st term | SW region | x-component)",
			"Lamb vect. circ. mom. coeff. VOLUME form. (1st term | VISCOUS region | x-component)",
			"Lamb vect. circ. mom. coeff. VOLUME form. (1st term | z-component)",
			"Lamb vect. circ. mom. coeff. VOLUME form. (2nd term | SW region | x-component)",
			"Lamb vect. circ. mom. coeff. VOLUME form. (2nd term | VISCOUS region | x-component)",
			"Lamb vect. circ. mom. coeff. VOLUME form. (2nd term | z-component)",
			"Liu Compress. corr. coeff. (1st term | x-component)",
			"Liu Compress. corr. coeff. (1st term | z-component)",
			"ONERA Compress. corr. coeff. (VOLUME formul. | WAVE region | x-component)",
			"ONERA Compress. corr. coeff. (VOLUME formul. | VISCOUS region | x-component)",
			"ONERA Compress. corr. coeff. (VOLUME formul. | SPURIOUS region | x-component)",
			"ONERA Compress. corr. coeff. (VOLUME formul. | z-component)",
			"Drag coeff. (TD method 1 - WAVE)",
			"Drag coeff. (TD method 1 - VISCOUS)",
			"Drag coeff. (TD method 1 - SPURIOUS)",
			"Entropy drag coeff. (TD method 1 - WAVE)",
			"Entropy drag coeff. (TD method 1 - VISCOUS)",
			"Entropy drag coeff. (TD method 1 - SPURIOUS)",
			"Lamb vect. circ. mom. coeff. (WAVE region surface integral)",
			"Lamb vect. circ. mom. coeff. (1st term | SPURIOUS region | x-component)",
			"Lamb vect. circ. mom. coeff. (2nd term | SPURIOUS region | x-component)",
			"Unsteady term coeff. (volume integral | x-component)",
			"Unsteady term coeff. (volume integral | z-component)",
			"Unsteady term coeff. (C-formulation volume integral | x-component)",
			"Unsteady term coeff. (C-formulation volume integral | z-component)",
			"Non-Inertial term coeff. (volume integral | x-component)",
			"Non-Inertial term coeff. (volume integral | z-component)",
			"Lamb vect. circ. mom. coeff. SURFACE form. (x-component)",
			"Lamb vect. circ. mom. coeff. SURFACE form. (z-component)",
			"Explicit viscous stress contribution coeff. (1st term | x-component)",
			"Explicit viscous stress contribution coeff. (1st term | z-component)",
			"Explicit viscous stress contribution coeff. (2nd term | x-component)",
			"Explicit viscous stress contribution coeff. (2nd term | z-component)",
			"ONERA compressibility correction coeff. (x-component)",
			"ONERA compressibility correction coeff. (z-component)",
			"Vortex force coeff. coeff. on body boundary (eff. inviscid flow | x-component )",
			"Vortex force coeff. coeff. on body boundary (eff. inviscid flow | z-component)",
			"Compress. corr. coeff. on body boundary (eff. inviscid flow | x-component)",
			"Compress. corr. coeff. on body boundary (eff. inviscid flow | z-component)",
			"Lamb vect. circ. mom. coeff. on body boundary (eff. inviscid flow | x-component)",
			"Lamb vect. circ. mom. coeff. on body boundary (eff. inviscid flow | z-component)",
			"Liu Compress. corr. coeff. (2nd term | x-component)",
			"Liu Compress. corr. coeff. (2nd term | z-component)",
			"Formulation C. coeff. (1st term | x-component)",
			"Formulation C. coeff. (1st term | z-component)",
			"Formulation C. coeff. (2nd term | x-component)",
			"Formulation C. coeff. (2nd term | z-component)",
			"Formulation C. coeff. (3rd term | x-component)",
			"Formulation C. coeff. (3rd term | z-component)",
			"Lamb vect. circ. mom. coeff. (VISCOUS region surface integral | x-component)",
			"Lamb vect. circ. mom. coeff. (VISCOUS region surface integral | z-component)",
			"ONERA compressibility corr. coeff. (VISCOUS region surface integral | x-component)",
			"ONERA compressibility corr. coeff. (VISCOUS region surface integral | z-component)",
			"Formulation C. coeff. (VISCOUS region surface integral | x-component)",
			"Formulation C. coeff. (VISCOUS region surface integral | z-component)",
			"Lamb vect. circ. mom. coeff. (VISCOUS WAKE surface integral | x-component)",
			"Lamb vect. circ. mom. coeff. (VISCOUS WAKE surface integral | z-component)",
			"ONERA compressibility corr. coeff. (VISCOUS WAKE surface integral | x-component)",
			"ONERA compressibility corr. coeff. (VISCOUS WAKE surface integral | z-component)",
			"Unsteady term coeff. on body surface (x-component)",
			"Unsteady term coeff. on body surface (z-component)",
			"Unsteady term coeff. on body surface (C-formulation | x-component)",
			"Unsteady term coeff. on body surface (C-formulation | z-component)",
			"Non-Inertial term coeff. (S_far integral | x-component)",
			"Non-Inertial term coeff. (S_far integral | z-component)",
			"Lamb vect. circ. mom. coeff. (WAVE region surface integral | x-component)",
			"Lamb vect. circ. mom. coeff. (WAVE region surface integral | z-component)",
			"ONERA compressibility corr. coeff. (WAVE region surface integral | x-component)",
			"ONERA compressibility corr. coeff. (WAVE region surface integral | z-component)",
			"Formulation C. coeff. (WAVE region surface integral | x-component)",
			"Formulation C. coeff. (WAVE region surface integral | z-component)"
		};
		string DBG_var_names[numDBG];
		vector < vector < vector < double > > > values_DBG(Nt, vector < vector < double > > (numDBG-1, vector < double > ((i_TOT - 1) * (kMax[0] - 1),0)));
		if (!debug_mode) {
			numDBG=0;
		}		
		if (debug_mode) {
			DBG_var_names[0] = "flag_ORG";
			DBG_var_names[1] = "Cd_P_vrt";
			DBG_var_names[2] = "Cd_irr_thd";
			DBG_var_names[3] = "Cd_i";
		}
		if (TEC_OUT==1) {
			cout << "WARNING:   Tecplot ASCII not implemented in unsteady version, sorry.\n";
		} else if (TEC_OUT == 2) {
			//Write Tecplot Binary Load-on-Demand file to output solution and post-processed data AS A SINGLE BLOCK (and ZONE).
			//Different time steps are attributed to different zones.
			void* outputFileHandle = NULL;
			int32_t outputDebugInfo = 0;
			if (debug_mode) {
				int32_t outputDebugInfo = 1;				
			}
			ostringstream outputStream;
				for (int32_t var = 0; var < numVars; var++) {
					outputStream << in_var_names[var];
					outputStream << ',';
				}
				for (int32_t var = 0; var < compVarsN + compVarsC + numAF + numFlags + numAF + numAF_ADD + numDBG; var++) {
					if (var < compVarsN + compVarsC + numAF + numFlags) {
						outputStream << out_var_names[var];
						outputStream << ',';
					} else if ( ((m_SEL==0)||(k_SEL==0)) && (var < compVarsN + compVarsC + numAF + numFlags + numAF + numAF_ADD) ) {
						outputStream << AF_var_names[var-(compVarsN + compVarsC + numAF + numFlags)];
						outputStream << ',';
					} else if (var >= compVarsN + compVarsC + numAF + numFlags + numAF + numAF_ADD) {
						outputStream << DBG_var_names[var - (compVarsN + compVarsC + numAF + numFlags + numAF + numAF_ADD)];
						if (var < compVarsN + compVarsC + numAF + numFlags + numAF + numAF_ADD + numDBG - 1) {
							outputStream << ',';
						}
					}
				}
			int32_t fileFormat = 1; // .szplt
			int numAF_OUT=numAF; int numAF_ADD_OUT=numAF_ADD;
				if ((m_SEL>0)&&(k_SEL>0)) {
					numAF_OUT=0;
					numAF_ADD_OUT=0;
				}
			int32_t varTypes[numVars+ compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG];
				for (int i = 0; i < numVars + compVarsN + compVarsC + numAF; i++) {
					varTypes[i] = 2;
				}
				for (int i = numVars + compVarsN + compVarsC + numAF; i < numVars + compVarsN + compVarsC + numAF + numFlags; i++) {
					varTypes[i] = 5;
				}
				for (int i = numVars + compVarsN + compVarsC + numAF + numFlags; i < numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT; i++) {
					varTypes[i] = 2;
				}
				for (int i = numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT; i < numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG; i++) {
					if (i==numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT) {
						varTypes[i] = 5;
					} else {
						varTypes[i] = 2;
					}
				}
			int32_t shareVarFromZone[numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG];
				for (int i = 0; i < numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG; i++) {
					shareVarFromZone[i] = 0;
				}
			int32_t valueLocation[numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG];
				for (int i = 0; i < numVars + compVarsN; i++) {
					valueLocation[i] = 1;
				}
				for (int i = numVars + compVarsN; i < numVars + compVarsN+ compVarsC + numAF + numFlags; i++) {
					valueLocation[i] = 0;
				}
				for (int i = numVars + compVarsN + compVarsC + numAF + numFlags; i < numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT; i++) {
					valueLocation[i] = 1;
				}
				for (int i = numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT; i < numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG; i++) {
					valueLocation[i] = 0;
				}
			int32_t passiveVarList[numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG];
				for (int i = 0; i < numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT + numAF_ADD_OUT + numDBG; i++) {
					passiveVarList[i] = 0;
				}
			int32_t shareConnectivityFromZone = 0;
			int64_t numFaceConnections = 0;
			int32_t faceNeighborMode = 0;
			int32_t outputZone;
			string TecFileName=argv[2];
				TecFileName.resize(TecFileName.length()-3);
				TecFileName=TecFileName+"szplt";
			R = tecFileWriterOpen(TecFileName.c_str(), "BF2D_DataSet", outputStream.str().c_str(), fileFormat, 0, 2, NULL, &outputFileHandle);
			R = tecFileSetDiagnosticsLevel(outputFileHandle, outputDebugInfo);
			vector < vector < vector < double > > > values_N(Nt, vector < vector < double > > (numVars + compVarsN, vector < double > (i_TOT * kMax[0],0)));
				for (int j = 0; j < numVars + compVarsN; j++) {
					for (int k = 0; k < kMax[0]; k++) {
						for (int m = 0; m < i_TOT; m++) {
							for (int t=0; t<Nt; t++) {
								values_N[t][j][i_TOT * k + m] = nodal_values[m][k][j][t];
							}
						}
					}
				}
				nodal_values.clear();
			vector < vector < vector < double > > > values_C(Nt, vector < vector < double > > (compVarsC, vector < double > ((i_TOT-1) * (kMax[0]-1),0)));
				for (int j = 0; j < compVarsC; j++) {
					for (int k = 0; k < kMax[0]-1; k++) {
						for (int m = 0; m < i_TOT-1; m++) {
							for (int t=0; t<Nt; t++) {
								values_C[t][j][(i_TOT-1) * k + m] = cell_values[m][k][numVars+compVarsN+j][t];
							}
						}
					}
				}
				cell_values.clear();
			vector < vector < vector < double > > > values_AF_elem(Nt, vector < vector < double > > (numAF, vector < double > ((i_TOT-1) * (kMax[0]-1),0)));
				for (int j = 0; j < numAF; j++) {
					for (int k = 0; k < kMax[0]-1; k++) {
						for (int m = 0; m < i_TOT-1; m++) {
							for (int t=0; t<Nt; t++) {
								values_AF_elem[t][j][(i_TOT-1) * k + m] = AF_elem[m][k][j][t];
							}
						}
					}
				}
			vector < vector < vector < uint8_t > > > values_F(Nt, vector < vector < uint8_t > > (numFlags, vector < uint8_t > ((i_TOT-1) * (kMax[0]-1),0)));
				for (int j = 0; j < numFlags; j++) {
					for (int k = 0; k < kMax[0]-1; k++) {
						for (int m = 0; m < i_TOT-1; m++) {
							for (int t=0; t<Nt; t++) {
								values_F[t][j][(i_TOT-1) * k + m] = flags[m][k][j][t];
							}
						}
					}
				}
				flags.clear();
			vector < vector < vector < double > > > values_AF;
				if ((m_SEL==0)||(k_SEL==0)) {
					values_AF = vector < vector < vector < double > > > (Nt, vector < vector < double > > (numAF + numAF_ADD, vector < double > (i_TOT * kMax[0],0)));
					for (int j = 0; j < numAF + numAF_ADD; j++) {
						for (int k = 0; k < kMax[0]; k++) {
							for (int m = 0; m < i_TOT; m++) {
								if ((k==kMax[0]-1)&&(m>=floor((i_TOT-1)/2)+1)&&(j==0)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Clp_nf[t]+Clf_nf[t];		//Set 1st aerodynamic force component at near field lift coeff. value at grid far-field boundary;
									}
								} else if ((k==kMax[0]-1)&&(m>=floor((i_TOT-1)/2)+1)&&(j==1)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Cdp_nf[t]+Cdf_nf[t];		//Set 2nd aerodynamic force component at near field drag coeff. value at grid far-field boundary;
									}
								} else if ((k==kMax[0]-1)&&(m>=floor((i_TOT-1)/2)+1)&&(j==2)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Clf_nf[t];        			//Set 3rd aerodynamic force component at near field friction lift coeff. value at grid far-field boundary;
									}
								} else if ((k==kMax[0]-1)&&(m>=floor((i_TOT-1)/2)+1)&&(j==3)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Cdf_nf[t];        			//Set 4th aerodynamic force component at near field friction drag coeff. value at grid far-field boundary;
									}
								} else if (m<floor((i_TOT-1)/2)+1) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = AF[m][k][j][t];
									}
								} else if ((m==i_TOT-1)&&(j==0)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Clp_nf[t]+Clf_nf[t];		//Set 1st aerodynamic force component at near field lift coeff. value at grid far-field boundary;
									}
								} else if ((m==i_TOT-1)&&(j==1)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Cdp_nf[t]+Cdf_nf[t];		//Set 2nd aerodynamic force component at near field drag coeff. value at grid far-field boundary;
									}
								} else if ((m==i_TOT-1)&&(j==2)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Cdf_nf[t];        			//Set 3rd aerodynamic force component at near field friction lift coeff. value at grid far-field boundary;
									}
								} else if ((m==i_TOT-1)&&(j==3)) {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = Cdf_nf[t];       			 //Set 4th aerodynamic force component at near field friction drag coeff. value at grid far-field boundary;
									}
								} else {
									for (int t=0; t<Nt; t++) {
										values_AF[t][j][i_TOT * k + m] = 0;
									}
								}
							}
						}
					}
				}
				AF.clear();
			vector < vector < uint8_t > > values_F_DBG(Nt, vector < uint8_t > ((i_TOT - 1) * (kMax[0] - 1), 0));
				for (int k = 0; k < kMax[0]-1; k++) {
					for (int m = 0; m < i_TOT-1; m++) {
						for (int t=0; t<Nt; t++) {
							values_F_DBG[t][(i_TOT-1) * k + m] = flag_org[m][k][t];
						}
					}
				}
				flag_org.clear();
				if (debug_mode) {
					for (int k = 0; k < kMax[0]-1; k++) {
						for (int m = 0; m < i_TOT-1; m++) {
							for (int t=0; t<Nt; t++) {
								values_DBG[t][0][(i_TOT-1) * k + m] = AF_elem[m][k][14][t] + AF_elem[m][k][15][t] + AF_elem[m][k][17][t] + AF_elem[m][k][18][t] + AF_elem[m][k][33][t] +
																	AF_elem[m][k][34][t] + AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t];
								values_DBG[t][1][(i_TOT-1) * k + m] = AF_elem[m][k][26][t] + AF_elem[m][k][27][t] + AF_elem[m][k][28][t];
								values_DBG[t][2][(i_TOT-1) * k + m] = AF_elem[m][k][0][t] + AF_elem[m][k][2][t] - (AF_elem[m][k][22][t] + AF_elem[m][k][23][t] + AF_elem[m][k][24][t]);
							}
						}
					}
				}
				AF_elem.clear();
			file=line.data();
			for (int t=0; t<Nt; t++) {
				stringstream middle;		
				middle << std::setw(pt-us-1) << std::setfill('0') << atof(line.substr(us+1, pt-us-1).data())+t;
				file=line.substr(0, us+1)+middle.str()+line.substr(pt);
				R = tecZoneCreateIJK(outputFileHandle, file.data(), i_TOT, 1, kMax[0], &varTypes[0],
					&shareVarFromZone[0], &valueLocation[0], &passiveVarList[0],
					shareConnectivityFromZone, numFaceConnections, faceNeighborMode, &outputZone);
				for (int j = 0; j < numVars + compVarsN; j++) {
					R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, j + 1, 0, i_TOT * kMax[0], &values_N[t][j][0]);
				}
				for (int j = 0; j < compVarsC; j++) {
					R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, numVars + compVarsN + j + 1, 0, (i_TOT - 1) * (kMax[0] - 1), &values_C[t][j][0]);
				}
				for (int j = 0; j < numAF; j++) {
					R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, numVars + compVarsN + compVarsC + j + 1, 0, (i_TOT - 1) * (kMax[0] - 1), &values_AF_elem[t][j][0]);
				}
				for (int j = 0; j < numFlags; j++) {
					R = tecZoneVarWriteUInt8Values(outputFileHandle, outputZone, numVars + compVarsN + compVarsC + numAF + j + 1, 0, (i_TOT - 1) * (kMax[0] - 1), &values_F[t][j][0]);
				}
				for (int j = 0; j < numAF_OUT+ numAF_ADD_OUT; j++) {
					R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, numVars + compVarsN + compVarsC + numAF + numFlags + j + 1, 0, i_TOT * kMax[0], &values_AF[t][j][0]);
				}
				for (int j = 0; j < numDBG; j++) {
					if (j==0) {
						R = tecZoneVarWriteUInt8Values(outputFileHandle, outputZone, numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT+ numAF_ADD_OUT + j + 1, 0, (i_TOT - 1) * (kMax[0] - 1), &values_F_DBG[t][0]);
					}else {			
						R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, numVars + compVarsN + compVarsC + numAF + numFlags + numAF_OUT+ numAF_ADD_OUT + j + 1, 0, (i_TOT - 1) * (kMax[0] - 1), &values_DBG[t][j-1][0]);
					}
				}
				R = tecZoneSetUnsteadyOptions(outputFileHandle, outputZone, time[t], t+1);
			}
			R = tecFileWriterClose(&outputFileHandle);
		}
	}
}
