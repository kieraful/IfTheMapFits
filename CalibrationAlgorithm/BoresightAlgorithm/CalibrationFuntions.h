#ifndef BORESIGHT_ITMF
#define BORESIGHT_ITMF

#include <iostream>
#include <Eigen\Dense> //for Eigen library
#include <Eigen\Core>
#include <fstream> //for read/write files
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iomanip> 
#include <string>
#include <cmath>  


//PCL includes
#include <pcl/features/3dsc.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/feature.h>
//#include <pcl/features/fpfh.h>
#include <pcl/features/impl/fpfh.hpp>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/io.h>
//#include <pcl/recognition/implicit_shape_model.h>
//#include <pcl/recognition/impl/implicit_shape_model.hpp>
#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/common/time.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/filters/statistical_outlier_removal.h>

using namespace std;
using namespace Eigen;

typedef Matrix<double, 2, 1> Matrix2b1;
typedef Matrix<double, Dynamic, 3> Matrixdb3;
typedef Matrix<double, 3, 3> Matrix3b3;
typedef Matrix<double, 3, 4> Matrix3b4;
typedef Matrix<double, Dynamic, 2> Matrixdby2;
typedef Matrix<int, Dynamic, 2> Matrixdby2i;

//typedef pcl::PointCloud<pcl::PointXYZI>::Ptr PointCloudXYZIptr;
//typedef pcl::PointCloud<pcl::PointXYZI> PointCloudXYZI;


typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloudXYZptr;
typedef pcl::PointCloud<pcl::PointXYZ> PointCloudXYZ;

#define MaxMatSize 10000000
#define PI 3.14159265358979323846 
//#define RAD2DEG = (180.0 / pi)


const double RAD2ARC = (3600 * 180 / PI);
const double RAD2DEG = (180 / PI);

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
	PointCloudXYZptr points_on_plane;

};

struct Scene {
	vector<Plane> planes;
	double X, Y, Z, omega, phi, kappa;
	

};

struct UniquePlanes {
	vector<RowVector3d> mapping_vec;
	vector<int> frequency;
	vector<Plane> unique_planes;

};

struct LidarPt {
	double x;
	double y;
	double z;
	double intensity;
	double time;
	
};

// ------------------------------ DR. SHAHBAHZI CODE ----------------------------


void Read_Mat(char *FileName, MatrixXd& m);//reads a formatted file to a matrix

void Write_Mat(char *FileName, MatrixXd & m, int decimal_precision); //writes a matrix to a formatted file

void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i); //Takes 3 Anges, makes 3x3 matrix

void Convert_R_to_Angles(Matrix3b3 R, double& Omega, double& Phi, double& Kappa); //Takes 3x3 matrix, makes 3 Anges

void  Normalization_Condition(MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& H1, MatrixXd& H2);

// ------------------------------------------------------------------------------

bool Read_Lidar_points(char * FileName, PointCloudXYZptr cloud);//reads a PCD file to a PCL point cloud

void Find_closest_points(int num_find_points, LidarPt point, MatrixXd search_points);

double euclidian_dist(double x1, double y1, double z1, double x2, double y2, double z2);

vector<Plane> FitPlanes(PointCloudXYZptr cloud_filtered, int max_planes = -1, bool make_files=false);

void filter_and_downsample(PointCloudXYZptr &input_cloud, float leaf_size = 0.001f);

void visualize_planes(vector<Plane> planes);

void visualize_cloud(PointCloudXYZptr cloud);

//void find_consensus_planes(vector<Plane> &all_planes);

bool sort_cloud(Plane plane_1, Plane plan_2);

bool sort_planes(Vector3d vec_1, Vector3d vec_2);

void save_planes(vector<Plane> planes);

void remove_outliers(PointCloudXYZptr &input_cloud, double search_n=50, double std_mul=1.0);

MatrixXd georeference_lidar_point(MatrixXd data, MatrixXd boresight_LA, MatrixXd boresight_angles);

UniquePlanes match_scenes(vector<Scene> scenes);

void remove_unfrequent(UniquePlanes &unique, int threshold=3);

void print_vector(vector<RowVector3d> print_vector);
void print_vector(vector<int> print_vector);

MatrixXd merge_data(MatrixXd IE_data, MatrixXd lidar_data);

#endif