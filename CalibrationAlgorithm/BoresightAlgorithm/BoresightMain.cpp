#include "CalibrationFuntions.h"


char * FILENAME = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\Data\\AppleWarehouse\\Orientation1.pcd";


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
	clog << "\n-------------------------STEP 1: Load PCD Scene Data-------------------------------------------------------\n";

	std::clog << "Opening file: " << FILENAME << " (can take up to 5 minutes)" << endl;
	PointCloudXYZptr Novatel_cloud(new PointCloudXYZ);
	if (!Read_Lidar_points(FILENAME, Novatel_cloud)) 
	{// Scene 1, Orientation 1
		clog << "\n\nProgram will close\n";
		return -1;
	}; 


	//clog << "\n-------------------------STEP 2: Mesh Data and Resample-------------------------------------------------------\n";

	// Done with cloud compare?

	clog << "\n-------------------------STEP 2: Filter Data-------------------------------------------------------\n";


	// Create the filtering object and downsample. (USE SUBSAMPLING INSTEAD)
	filter_and_downsample(Novatel_cloud, 0.1f);


	clog << "\n-------------------------STEP 3: Fit all planes-----------------------------------------------------\n";

	
	vector<Plane> planes_in_cloud = FitPlanes(Novatel_cloud);

	// Find the largest planes
	std::sort(planes_in_cloud.begin(), planes_in_cloud.end(), sort_cloud); // sort based off cloud size



	planes_in_cloud.resize(5); //truncate to keep largest planes
	//save planes
	save_planes(planes_in_cloud);
	

	clog << "\n-------------------------STEP 4: Downsample pts on Planes----------------------------------------------------\n";

		//TODO: downsample all points on each plane. These will be # of EQUATIONS
	clog << "Downsampling.....\n";

	for (int i = 0; i < planes_in_cloud.size(); i++) {
		remove_outliers(planes_in_cloud[i].points_on_plane, 100, 1);
		filter_and_downsample(planes_in_cloud[i].points_on_plane, 1.0f);
	}
	

	// --------------------------------------------------------------------------------------------------------------------------------




	// TODO : GET GNSS/INS integrated data (Integrated from Inertial Explorer)

	//Read in IE output and save to matrix
	MatrixXd GNSS_INS_data;
	Read_Mat("LiDAR_Georeferencing_Output_noheader.txt", GNSS_INS_data);

	// TODO: Match timestamps of GNSS+INS data to that of LiDAR data



	clog << "\n-------------------------STEP 4: Visualize-----------------------------------------------------\n";

	// VISUALIZE
	visualize_planes(planes_in_cloud);


	return 0;
}