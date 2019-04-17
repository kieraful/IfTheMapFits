#ifndef LS_H
#define LS_H
#endif

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <Eigen\Dense> //for Eigen library
#include <Eigen\Core>
#include <fstream> //for read/write files
#include <algorithm>
#include <stdexcept>
#include <string>
#include <iomanip> 

#define MaxMatSize 10000000
#define PI 3.14159265358979323846


using namespace std;
using namespace Eigen;


void Write_Mat(char *FileName, MatrixXd &m, int decimal_precision);

void Read_Mat(char *FileName, MatrixXd& m);

//Function to fill a matrix with info from vector of vectors
void vec2mat(vector<RowVectorXd>& vec, MatrixXd& mat, int cols);

//Function to compute A
void computeAandw(MatrixXd& A, MatrixXd& H, MatrixXd& w, MatrixXd& V, int u, int numPlanes, int numScans, int numLidPts, MatrixXd& bs_params, MatrixXd& obs_bs, MatrixXd& estscans, MatrixXd& plane_details, MatrixXd& scene_details, MatrixXd& point_details);

//Function to compute 6 rows of A for a scan
MatrixXd computeA1(int u, int numPlanes, int scanNum);

//Function to compute a row of A for a plane
MatrixXd computeA2(int u, int planeNum, double n_xpg, double n_ypg, double n_zpg);

//Function to compute elements of a rotation matrix
MatrixXd RotMatElements(double w, double phi, double K);

//Function to compute the derivatives of a point equation wrt rotation matrix elements for Rbjg
MatrixXd PtEqnWrtRotbjg(double x_Sjb, double y_Sjb, double z_Sjb, double w_Sb, double phi_Sb, double K_Sb, double x_sj, double y_sj, double z_sj, double n_xpg, double n_ypg, double n_zpg);

//Function to compute the derivatives of a point equation wrt rotation matrix elements for RSb
MatrixXd PtEqnWrtRotSb(double w_bjg, double phi_bjg, double K_bjg, double x_sj, double y_sj, double z_sj, double n_xpg, double n_ypg, double n_zpg);

//Function to compute derivatives of rotation matrix elements wrt rotation angles
MatrixXd RotWrtAngles(double w, double phi, double K);

//Function to compute a row of A for a lidar point
MatrixXd computeA3(int u, int numPlanes, int planeNum, int scanNum,
	double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,	
	double n_xpg, double n_ypg, double n_zpg,	
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj);

//Function to compute 6 rows of w for a scan
MatrixXd computew1(double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_GPS, double y_GPS, double z_GPS, 
	double w_INS, double phi_INS, double K_INS);

//Function to compute a row of w for a plane
MatrixXd computew2(double n_xpg, double n_ypg, double n_zpg);

//Function to compute a row of w for a lidar point
MatrixXd computew3(double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,	
	double n_xpg, double n_ypg, double n_zpg, double d_p,	
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj);
	
MatrixXd computeB1(int n, int numLidPts, int scanNum);
	
MatrixXd computeB3(int n, int ptNum, double n_xpg, double n_ypg, double n_zpg, double w_bjg, double phi_bjg, double K_bjg, double w_Sb, double phi_Sb, double K_Sb);
