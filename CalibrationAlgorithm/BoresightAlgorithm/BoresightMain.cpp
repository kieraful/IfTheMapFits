#include "CalibrationFuntions.h"




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

	//Initialize variables 
	vector<char *> pcd_files;
	char * file1 = "..\\..\\..\\Data\\AppleWarehouse\\Orientation1.pcd";
	pcd_files.push_back(file1); //vector of input files

	vector<Plane> planes;
	vector<Scene> scenes;

	for (int i = 0; i < pcd_files.size(); i++)
	{
		Scene temp_scene;

		// ---------------------------------------STEP 1: Load PCD Scene Data-----------------------------------------------------------------------------------------
		clog << "\n-------------------------STEP 1: Load PCD Scene Data-------------------------------------------------------\n";

		std::clog << "Opening file: " << pcd_files[i] << " (can take up to 5 minutes)" << endl;
		PointCloudXYZptr Novatel_cloud(new PointCloudXYZ);
		if (!Read_Lidar_points(pcd_files[i], Novatel_cloud))
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


		planes = FitPlanes(Novatel_cloud);



		// Find the largest planes
		std::sort(planes.begin(), planes.end(), sort_cloud); // sort based off cloud size



		planes.resize(5); //truncate to keep largest planes
		//save planes
		save_planes(planes);


		clog << "\n-------------------------STEP 4: Downsample pts on Planes----------------------------------------------------\n";

		//TODO: downsample all points on each plane. These will be # of EQUATIONS
		clog << "Downsampling.....\n";

		for (int i = 0; i < planes.size(); i++) {
			remove_outliers(planes[i].points_on_plane, 100, 1);
			filter_and_downsample(planes[i].points_on_plane, 1.0f);
		}


		// --------------------------------------------------------------------------------------------------------------------------------




		// TODO : GET GNSS/INS integrated data (Integrated from Inertial Explorer)

		//GNSS
		temp_scene.X = 0;
		temp_scene.Y = 0;
		temp_scene.Z = 0;
		//INS
		temp_scene.omega = 0;
		temp_scene.phi = 0;
		temp_scene.kappa = 0;
		
		// VISUALIZE
		visualize_planes(planes);

		//Add to scene
		temp_scene.planes = planes;
		
		scenes.push_back(temp_scene);
	}

	return 0;
}