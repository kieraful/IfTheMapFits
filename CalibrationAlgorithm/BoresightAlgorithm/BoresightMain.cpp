#include "CalibrationFuntions.h"

char * FILENAME = "..\\..\\..\\Data\\FirstDataset\\All_points_resampled.pcd";

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
	PointCloudXYZptr Novatel_cloud(new PointCloudXYZ);
	Read_Lidar_points(FILENAME, Novatel_cloud); // Scene 1, Orientation 1


	clog << "\n-------------------------STEP 2: Mesh Data and Resample-------------------------------------------------------\n";

	// Done with cloud compare?

	clog << "\n-------------------------STEP 2: Filter Data-------------------------------------------------------\n";

		//TODO: use PCL to filter data

	// Create the filtering object and downsample.
	PointCloudXYZptr filter_cloud = filter_and_downsample(Novatel_cloud, 0.1f);


	clog << "\n-------------------------STEP 3: Fit all planes-----------------------------------------------------\n";


		//TODO: Incorporate Plane fitting algorithm.
	
	vector<Plane> planes_in_cloud = FitPlanes(filter_cloud, 3, true);


		//TODO: Find how to uniquely describe planes, as output from plane-fitting

	
	// ---------------------------------------STEP 4: Downsample pts on Planes-----------------------------------------------------------------------------------------


		//TODO: downsample all points on each plane. These will be # of EQUATIONS

	// --------------------------------------------------------------------------------------------------------------------------------




	// TODO : GET GNSS/INS integrated data (Integrated from Inertial Explorer)





	// VISUALIZE
	visualize_planes(planes_in_cloud);


	return 0;
}