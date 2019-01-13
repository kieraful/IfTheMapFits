#ifndef BORESIGHT_ITMF
#define BORESIGHT_ITMF

#include <iostream>
#include <Eigen\Dense> //for Eigen library
#include <Eigen\Core>
#include <fstream> //for read/write files
#include <stdlib.h>// for system("pause") and exit_failure
#include <vector>
#include <math.h>  
#include <algorithm>
#include <stdexcept>
#include <iomanip> 
#include <string>

//PCL includes
//#include <pcl/io/pcd_io.h>
//#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>


using namespace std;
using namespace Eigen;

typedef Matrix<double, 2, 1> Matrix2b1;
typedef Matrix<double, Dynamic, 3> Matrixdb3;
typedef Matrix<double, 3, 3> Matrix3b3;
typedef Matrix<double, 3, 4> Matrix3b4;
typedef Matrix<double, Dynamic, 2> Matrixdby2;
typedef Matrix<int, Dynamic, 2> Matrixdby2i;


#define MaxMatSize 10000000
#define pi 3.14159265358979323846 

struct CameraParam {
	double PS; //pixel size
	double f_l; //focal length
	double xpp; //principal point in x
	double ypp; //principal point in y
	double K1; //Distortion parameter
	double K2;//Distortion parameter
	double K3;//Distortion parameter
	double P1;//Distortion parameter
	double P2;//Distortion parameter
	double S1;//Distortion parameter
	double S2;//Distortion parameter
	double Cn;
	double Rn;
	double sigma_obs; //observation standard deviation

};

struct Plane {

	double a1, a2, a3, b; // plane parameters
	MatrixXd points_in;

};

struct LidarPt {
	double x;
	double y;
	double z;
	double intensity;
	double time;
	
};

// ------------------------------ DR. SHAHBAHZI CODE ----------------------------

void Read_Lidar_points(char *FileName, MatrixXd& m);//reads a formatted file to a matrix

void Read_Mat(char *FileName, MatrixXd& m);//reads a formatted file to a matrix

void Write_Mat(char *FileName, MatrixXd & m, int decimal_precision); //writes a matrix to a formatted file

void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i); //Takes 3 Anges, makes 3x3 matrix

void Convert_R_to_Angles(Matrix3b3 R, double& Omega, double& Phi, double& Kappa); //Takes 3x3 matrix, makes 3 Anges

void  Normalization_Condition(MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& H1, MatrixXd& H2);

void Find_closest_points(int num_find_points, LidarPt point, MatrixXd search_points);

double euclidian_dist(double x1, double y1, double z1, double x2, double y2, double z2);


// ------------------------------------------------------------------------------






#endif