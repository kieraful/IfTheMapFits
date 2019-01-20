#include "CalibrationFuntions.h"

char * FILENAME = "..\\..\\..\\Data\\FirstDataset\\All_points_dec5.pcd";

int main() {
	/*
	--------------------------- IF THE MAP FITS BORESIGHT CALIBRATION - MAIN  ------------------------------------------------------

	

	--------------------------------------------------------------------------------------------------------------------------------

	*/
	// Program Header
	fprintf(stdout, "\n----------------------IF THE MAP FITS BORESIGHT CALIBRATION----------------------------\n");
	fprintf(stdout, "\n\tThe purpose of this program is to compute the Boresight calibration parameters of \n");
	fprintf(stdout, "\ta system with GNSS/INS/LiDAR sensors \n\n\n");
	
	fprintf(stdout, "    ___________________         \n");
	fprintf(stdout, "   |,-----.,-----.,---.\        \n");
	fprintf(stdout, "   ||     ||     ||    \\       \n");
	fprintf(stdout, "   |`-----'|-----||-----\`----. \n");
	fprintf(stdout, "   [       |    -||-   _|    (| \n");
	fprintf(stdout, "   [  ,--. |_____||___/.--.   | \n");
	fprintf(stdout, "   =-(( `))-----------(( `))-== \n");
	fprintf(stdout, "      `--'             `--'     \n");


	fprintf(stdout, "\n---------------------------------------------------------------------------------------\n");

	// ---------------------------------------STEP 1: Load PCD Scene Data-----------------------------------------------------------------------------------------

	std::clog << "Opening file: " << FILENAME << " (can take up to 5 minutes)" << endl;
	PointCloudXYZIptr Novatel_cloud(new PointCloudXYZIptr);
	Read_Lidar_points(FILENAME, Novatel_cloud); // Scene 1, Orientation 1


	// ---------------------------------------STEP 2: Median Filter Data-----------------------------------------------------------------------------------------

		//TODO: use PCL to filter data

	// Create the filtering object and downsample.
	pcl::VoxelGrid<pcl::PointXYZI> vox_grid;
	vox_grid.setInputCloud(Novatel_cloud);
	vox_grid.setLeafSize(0.01f, 0.01f, 0.01f);
	vox_grid.filter(*Novatel_cloud);


	// ---------------------------------------STEP 3: Fit all planes-----------------------------------------------------------------------------------------


		//TODO: Incorporate Plane fitting algorithm.
		//TODO: Find how to uniquely describe planes, as output from plane-fitting

	
	// ---------------------------------------STEP 4: Downsample pts on Planes-----------------------------------------------------------------------------------------


		//TODO: downsample all points on each plane. These will be # of EQUATIONS

	// --------------------------------------------------------------------------------------------------------------------------------




	// TODO : GET GNSS/INS integrated data (Integrated from Inertial Explorer)





	// VISUALIZE

	visualize_cloud(Novatel_cloud);







		return 0;
	}