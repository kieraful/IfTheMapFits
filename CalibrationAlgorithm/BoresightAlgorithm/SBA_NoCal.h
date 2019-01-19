#ifndef SBA_NoCal_H
#define SBA_NoCal_H

#include "CalibrationFuntions.h"

using namespace std;
using namespace Eigen;

#define MaxMatSize 10000000
#define pi 3.14159265358979323846 
#define EIGEN_USE_LAPACKE 
#define EIGEN_USE_LAPACKE_STRICT	


typedef Matrix<double, 2, 1> Matrix2b1;
typedef Matrix<double, Dynamic, 3> Matrixdb3;
typedef Matrix<double, 3, 3> Matrix3b3;
typedef Matrix<double, 3, 4> Matrix3b4;
typedef Matrix<double, Dynamic, 2> Matrixdby2;
typedef Matrix<int, Dynamic, 2> Matrixdby2i;
typedef Matrix<bool, Dynamic, Dynamic> MatrixBools;
// Read and write function are modified versions of DMP_BBO library
void Read_Mat(char *FileName, MatrixXd& m);//reads a formatted file to a matrix
void Write_Mat(char *FileName, MatrixXd &m, int decimal_precision); //writes a matrix to a formatted file

																	//3x3 rotation matrix defined as Rz(Kappa)*Ry(Phi)*Rx(Omega), rotates the image space to become parallel to the object space
void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i);

//cehcks whether number k belongs to vector V and returns the index I that shows where in vector V the number K exists.
bool does_exist(vector<double> V, double k, int & I);
bool does_exist_eigen(VectorXi V, int k, int & I);


void removeRow(MatrixXd& matrix, unsigned int rowToRemove);//remove a row from a matrix; rowToRemove is the row number
void removeColumn(MatrixXd& matrix, unsigned int colToRemove);//remove a column from a matrix; colToRemove is the column number


															  // to sort a matrix rows based on a given column number 
															  // Simple little templated comparison function 
															  //source stackoverflow.com/questions/39693909/sort-eigen-matrix-column-values-by-ascending-order-of-column-1-values
template <typename MatrixT>
bool compareRows(MatrixT a, MatrixT b) {
	return a(0, 0) < b(0, 0);
}
// These are the 6 template arguments to every Eigen matrix
template <typename Scalar, int rows, int cols, int options, int maxRows, int maxCols>
void  SortRows(Matrix<Scalar, rows, cols, options, maxRows, maxCols> & sorted, Matrix<Scalar, rows, cols, options, maxRows, maxCols>& target, int coltofocus) {

	// Manually construct a vector of correctly-typed matrix rows
	Matrix<Scalar, rows, cols, options, maxRows, maxCols> target_temp;
	target_temp = target;

	target_temp.col(0) = target.col(coltofocus);
	target_temp.col(coltofocus) = target.col(0);

	vector<Eigen::Matrix<Scalar, 1, cols>> matrixRows;
	for (unsigned int i = 0; i < target_temp.rows(); i++) {
		matrixRows.push_back(target_temp.row(i));
	}


	sort(matrixRows.begin(),
		matrixRows.end(),
		compareRows<Matrix<Scalar, 1, cols>>);


	Matrix<Scalar, rows, cols, options, maxRows, maxCols> sorted_temp;
	sorted_temp = target;

	for (unsigned int i = 0; i < matrixRows.size(); i++) {
		sorted_temp.row(i) << matrixRows[i];
	}

	sorted = sorted_temp;
	sorted.col(0) = sorted_temp.col(coltofocus);
	sorted.col(coltofocus) = sorted_temp.col(0);

	return;
};

// a base simple class to save the camera EOPs and IOPs, update them, write them to file.
class CameraParam
{
private:
	//none
public:

	// calibration (IO) parameters for model 0
	MatrixXd  xpp_0; //principal point x-offset (mm)
	MatrixXd  ypp_0; //principal point y-offset (mm)
	MatrixXd  f_l_0; //principal distance (focal length) (mm)
	MatrixXd  K1_0;
	MatrixXd  K2_0;
	MatrixXd  K3_0;
	MatrixXd  P1_0;
	MatrixXd  P2_0;
	MatrixXd  S1_0;
	MatrixXd  S2_0;
	MatrixXd  PS_0; //pixel size  (mm)
	MatrixXd  Cn_0; //number of columns of the image (pixels)
	MatrixXd  Rn_0; //number of rows of the image (pixels)

					// calibration (IO) parameters for model 1
	MatrixXd  cx_1;
	MatrixXd  cy_1;
	MatrixXd  fx_1;
	MatrixXd  fy_1;
	MatrixXd  K1_1;
	MatrixXd  K2_1;
	MatrixXd  K3_1;
	MatrixXd  P1_1;
	MatrixXd  P2_1;

	// calibration (IO) parameters for model 2
	MatrixXd  cx_2;
	MatrixXd  cy_2;
	MatrixXd  C_2;
	MatrixXd  D_2;
	MatrixXd  E_2;
	MatrixXd  F_2;
	MatrixXd  P2_2;
	MatrixXd  P3_2;
	MatrixXd  P4_2;

	VectorXi Camera_labels;
	MatrixBools Calibrated_0, Calibrated_1, Calibrated_2;

	int num_sensors_model0; //number of sensors that follow modelid 0
	int num_sensors_model1; //number of sensors that follow modelid 1
	int num_sensors_model2; //number of sensors that follow modelid 2

	MatrixXd sensor_ID; //each image is taken by which sensor ; numbered from zero in each model.
	MatrixXd model_ID; //each image is following which model: 
					   //if traditional pinhole model, then modelid=0
					   // if OpenCV pinhole model, then modelid=1
					   // if equi-distance Pix4D fisheye model, then modelid=2;
					   //for example if you have 6 sensors and one image is taken by each, 1 one of them follows pinhole, 1 of them opencv, 4 of them fisheye
					   //then you sensor_IDs are 0, 1, 2, 3 , 4, 5
					   //you model_IDs are 0, 1, 2, 2, 2, 2
	Matrixdb3 Bundestim_ori; //EO parameters format: Rmatrix;Tvector;Omega Phi Kappa

	CameraParam(int nsensorwithmodelid0, int nsensorwithmodelid1, int nsensorwithmodelid2); // constructor
	~CameraParam(); // destructor 

	void Re_init_Estimori(int i, double Omega, double Phi, double Kappa, double Xo, double Yo, double Zo); //updates EO parameters
	void ResetIOPs2(int i, bool calibrated, double P2v, double P3v, double P4v, double cxv, double cyv, double Cv, double Dv, double Fv);
	void ResetIOPs1(int i, bool calibrated, double K1v, double K2v, double K3v, double P1v, double P2v, double cxv, double cyv, double fxv, double fyv);
	void ResetIOPs0(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double P1v, double P2v, double S1v, double S2v, double xppv, double yppv, double f_lv);
	void ResetEOPs(Matrixdb3& Estimori); //updates EO and IO parameters of a camera object
	void WriteCameraToFile(char *FileName, char *FileName2); //writes all camera parameters to a file

};


void Data_Preparation(char* FeatSavePP, char* DirectEO,
	char* ALLTIEPOINTS, char* AllModels,
	VectorXi & Point_Labels_Final,
	Matrixdby2 & Obspixels_p, MatrixXd & Sigmaxy_p, Matrixdby2i & Index_p_i,
	VectorXi & Adress_i_p, VectorXi & Adress_p_i,
	VectorXi & Index_i_p, CameraParam& camera_params); //read the data from text files and save them to the sparsely formatted matrices (see the supplementary notes)

void Data_Preparation_Camonly(char* FeatSavePP, char* DirectEO, char* AllModels, CameraParam& camera_params);


//perform the SBA!
void Perform_BA(VectorXi & Point_Labels,
	Matrixdby2 & Obspixels_p, MatrixXd & Sigmaxy_p, Matrixdby2i & Index_p_i,
	VectorXi & Adress_i_p, VectorXi & Adress_p_i,
	VectorXi & Index_i_p, CameraParam& camera_params, MatrixXd& ground, MatrixXd delta_trans, double IO_scalfact, MatrixXd& Sigma_Xcap_C, MatrixXd& Sigma_Xcap_G);



#endif