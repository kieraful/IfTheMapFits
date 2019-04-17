#ifndef BAFUNCTIONS_H
#define BAFUNCTIONS_H
#endif

#include "stdafx.h"
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

using namespace std;
using namespace Eigen;


//Function to fill a matrix with info from vector of vectors
void vec2mat(vector<RowVectorXd>& vec, MatrixXd& mat, int cols);

//Function to compute A
void computeAandw(MatrixXd& A_full, MatrixXd& A, MatrixXd& H, MatrixXd& w_full, MatrixXd& w, MatrixXd& V, int u, int numPlanes, int numScans, int numLidPts, MatrixXd bs_params, MatrixXd plane_details, MatrixXd scene_details, MatrixXd point_details, MatrixXd GNSS_INS_data);

//Function to compute 6 rows of A for a scan
MatrixXd computeAScan(int u, int numPlanes, int scanNum);

//Function to compute a row of A for a plane
MatrixXd computeAPlane(int u, int planeNum, double n_xpg, double n_ypg, double n_zpg);

//Function to compute elements of a rotation matrix
MatrixXd RotMatElements(double w, double phi, double K);

//Function to compute the derivatives of a point equation wrt rotation matrix elements for Rbjg
MatrixXd PtEqnWrtRotbjg(double x_Sjb, double y_Sjb, double z_Sjb, double w_Sb, double phi_Sb, double K_Sb, double x_sj, double y_sj, double z_sj, double n_xpg, double n_ypg, double n_zpg);

//Function to compute the derivatives of a point equation wrt rotation matrix elements for RSb
MatrixXd PtEqnWrtRotSb(double w_bjg, double phi_bjg, double K_bjg, double x_sj, double y_sj, double z_sj, double n_xpg, double n_ypg, double n_zpg);

//Function to compute derivatives of rotation matrix elements wrt rotation angles
MatrixXd RotWrtAngles(double w, double phi, double K);

//Function to compute a row of A for a lidar point
MatrixXd computeAPt(int u, int numPlanes, int planeNum, int scanNum,
	double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,	
	double n_xpg, double n_ypg, double n_zpg,	
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj);

//Function to compute 6 rows of w for a scan
MatrixXd computewScan(double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_GPS, double y_GPS, double z_GPS, 
	double w_INS, double phi_INS, double K_INS);

//Function to compute a row of w for a plane
MatrixXd computewPlane(double n_xpg, double n_ypg, double n_zpg);

//Function to compute a row of w for a lidar point
MatrixXd computewPt(double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,	
	double n_xpg, double n_ypg, double n_zpg, double dp_sum,	
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj);