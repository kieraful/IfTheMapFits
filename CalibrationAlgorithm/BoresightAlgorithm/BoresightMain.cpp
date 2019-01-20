#include "CalibrationFuntions.h"

char * FILENAME = "..\\..\\..\\Data\\FirstDataset\\All_points_dec5.pcd";

int main() {
	/*
	--------------------------- IF THE MAP FITS BORESIGHT CALIBRATION - MAIN  ------------------------------------------------------


	--------------------------------------------------------------------------------------------------------------------------------

	*/


	// ---------------------------------------STEP 1: Load PCD Scene Data-----------------------------------------------------------------------------------------

	pcl::PointCloud<pcl::PointXYZI> Novatel_cloud;
	Read_Lidar_points(FILENAME, Novatel_cloud); // Scene 1, Orientation 1

	// ---------------------------------------STEP 2: Median Filter Data-----------------------------------------------------------------------------------------

		//TODO: use PCL to filter data

	// ---------------------------------------STEP 3: Fit all planes-----------------------------------------------------------------------------------------


		//TODO: Incorporate Plane fitting algorithm.
		//TODO: Find how to uniquely describe planes, as output from plane-fitting

	
	// ---------------------------------------STEP 4: Downsample pts on Planes-----------------------------------------------------------------------------------------


		//TODO: downsample all points on each plane. These will be # of EQUATIONS

	// --------------------------------------------------------------------------------------------------------------------------------




	// TODO : GET GNSS/INS integrated data (Integrated from Inertial Explorer)





	// VISUALIZE

	//visualize_cloud(FILENAME);







		return 0;
	}