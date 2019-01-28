#include "CalibrationFuntions.h"

char * FILENAME = "..\\..\\..\\Data\\FirstDataset\\All_points_dec5.pcd";

char * FILENAME2 = "Upsampled_PointCloud.pcd";

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

	std::clog << "Opening file: " << FILENAME2 << " (can take up to 5 minutes)" << endl;
	PointCloudXYZptr Novatel_cloud2(new PointCloudXYZ);
	Read_Lidar_points(FILENAME2, Novatel_cloud2); // Scene 1, Orientation 1

	clog << "\n-------------------------STEP 2: Mesh Data-------------------------------------------------------\n";

	// Normal estimation*
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimation;
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud(Novatel_cloud);
	normal_estimation.setInputCloud(Novatel_cloud);
	normal_estimation.setSearchMethod(tree);
	normal_estimation.setKSearch(20);
	normal_estimation.compute(*normals);
	//* normals should not contain the point normals + surface curvatures

	clog << "\n-------------------------STEP 2: Filter Data-------------------------------------------------------\n";

		//TODO: use PCL to filter data

	// Create the filtering object and downsample.
	//PointCloudXYZIptr filter_cloud = filter_and_downsample(Novatel_cloud, 0.1f);


	clog << "\n-------------------------STEP 3: Fit all planes-----------------------------------------------------\n";


		//TODO: Incorporate Plane fitting algorithm.
	
	//vector<Plane> planes_in_cloud = FitPlanes(filter_cloud, 4);


		//TODO: Find how to uniquely describe planes, as output from plane-fitting

	
	// ---------------------------------------STEP 4: Downsample pts on Planes-----------------------------------------------------------------------------------------


		//TODO: downsample all points on each plane. These will be # of EQUATIONS

	// --------------------------------------------------------------------------------------------------------------------------------




	// TODO : GET GNSS/INS integrated data (Integrated from Inertial Explorer)





	// VISUALIZE
	//visualize_planes(planes_in_cloud);
	visualize_cloud(Novatel_cloud);
	visualize_cloud(Novatel_cloud2);


	return 0;
}