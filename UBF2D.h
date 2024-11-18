#include <vector>
#include "/usr/local/apps/tecplot2022r1/360ex_2022r1/include/TECIO.h"	//Required by TecIO

using namespace std;

void cell_metrics(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values, bool NS, int U_FORM, bool MOV_GRD, int t, int m, int k,
				  int numVars, int compVarsN, vector < bool > &near_ID, int i_TOT, int64_t* kMax, int grad_scheme, bool flux_scheme, bool ParForm, vector < double > &X_POLE, vector < double > &Z_POLE,
				  double& rx_1, double& rz_1, double& rx_2, double& rz_2, double& rx_3, double& rz_3, double& rx_4, double& rz_4,
				  double& dx_1, double& dz_1, double& dx_2, double& dz_2, double& dx_3, double& dz_3, double& dx_4, double& dz_4,
				  double& a_1, double& a_2, double& a_3, double& a_4, vector < vector < double > > &WLS);

void updateCompCorr(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values, int rho_id,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, int U_FORM, bool MOV_GRD, int mrho_form, bool flux_scheme, vector < bool > &near_ID,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem);

void updateFirstTerm(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values, int rho_id,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, int U_FORM, bool MOV_GRD, int mrho_form, bool flux_scheme, vector < bool > &near_ID,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem);

void updateFourthTerm(vector < vector < vector < vector < double > > > > &nodal_values,vector < vector < vector < vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < double > > > > &cell_values,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, bool flux_scheme, vector < bool > &near_ID,
					double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < bool > > > > &flags,
					vector < vector < vector < vector < double > > > > &AF_elem);

void updateParasite(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector <vector < double > > > > &cell_values,
					vector < vector < vector < vector < double > > > > &aux_values, vector < vector < vector <vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < bool > > > > &flags,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, int U_FORM, bool MOV_GRD, bool flux_scheme, bool smart_wds,
					vector < bool > &near_ID, double CdP_tsh, double rho_inf, double V_inf, double chord,
					double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem);

void updateParasiteNew(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector <vector < double > > > > &cell_values,
					vector < double > &X_POLE, vector < double > &Z_POLE, int rho_id, int vx_id, int vz_id, vector < vector < vector <vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < bool > > > > &flags,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, bool MOV_GRD, bool flux_scheme, bool smart_wds,
					vector < bool > &near_ID, double CdP_tsh, double rho_inf, double V_inf, double chord,
					double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem);

void updateF_grad_rho(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values,
					vector < vector < vector < vector < double > > > > &aux_values, vector < vector < vector < vector < double > > > > &aux_values_N,
					vector < vector < vector < vector < bool > > > > &flags,
					int t, int m, int k, int numVars, int compVarsN, int i_TOT, int64_t* kMax, bool NS, int U_FORM, bool MOV_GRD, bool flux_scheme, bool smart_wds,
					vector < bool > &near_ID, double CdP_tsh, double rho_inf, double V_inf,
					double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				  	double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				  	double a_1, double a_2, double a_3, double a_4, vector < vector < vector < vector < double > > > > &AF_elem);

void compute_AF(vector < vector < vector < vector < double > > > > &nodal_values, vector < vector < vector < vector < double > > > > &cell_values,
				vector < vector < vector < vector < double > > > > &aux_values, vector < vector < vector < vector < double > > > > &aux_values_N,
				vector < vector < vector < vector < bool > > > > &flags, vector < vector < vector < vector < double > > > > &AF_elem,
				int t, int mD, int kD, int numVars, int compVarsN, int numAF, int numAF_ADD, int i_TOT, int64_t* kMax, bool NS, bool flux_scheme, int grad_scheme, bool ParForm, bool limit_integr,
				vector < bool > &near_ID, double rho_inf, double V_inf, double chord, int rho_id, int vx_id, int vz_id, int vgx_id, int vgz_id, vector < double > &X_POLE, vector < double > &Z_POLE,
				int U_FORM, bool FT_FORM, double x_CoR, double z_CoR, double ampl_x, double freq_x, double t0_x, double ampl_z, double freq_z, double t0_z,
				double ampl_a, double freq_a, double t0_a, bool MOV_GRD, vector < double > &time,
				double rx_1, double rz_1, double rx_2, double rz_2, double rx_3, double rz_3, double rx_4, double rz_4,
				double dx_1, double dz_1, double dx_2, double dz_2, double dx_3, double dz_3, double dx_4, double dz_4,
				double a_1, double a_2, double a_3, double a_4, vector < vector < double > > &WLS, vector < vector < double > > &AERO);
				