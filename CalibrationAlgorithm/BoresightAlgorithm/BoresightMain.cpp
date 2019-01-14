#include "CalibrationFuntions.h"

char * FILENAME = "..\\..\\..\\Data\\FirstDataset\\All_points_dec5.pcd";

int main() {
	/*
	--------------------------- IF THE MAP FITS BORESIGHT CALIBRATION - MAIN  ------------------------------------------------------


	--------------------------------------------------------------------------------------------------------------------------------

	*/

	// TODO : GET LiDAR DATA (Either only feature data or full dataset, feature preferred) 
	pcl::PointCloud<pcl::PointXYZI> Novatel_cloud;
	Read_Lidar_points(FILENAME, Novatel_cloud);

	//Mesh Shapes
	//pcl::NormalEstimation<pcl::PointXYZI, pcl::Normal> normal_estimator;



	// TODO : GET GNSS/INS integrated data (Integrated from Inertial Explorer)





	// VISUALIZE

	visualize_cloud(FILENAME);







		return 0;
	}