#include "CalibrationFuntions.h"
//#include "BAFunctions.h"


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
char * file1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation1_b.pcd";
char * file2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation2_b.pcd";
char * file3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation3_b.pcd";
char * file4 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation4_b.pcd";

//char * file1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Orientation1_b.pcd";
//char * file2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Orientation2_b.pcd";
//char * file3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Orientation3_b.pcd";

//char * file1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
//char * file2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points_dec5.pcd";
//char * file3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";

//char * file1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
//char * file2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
//char * file3 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";

//pcd_files.push_back(file1); //vector of input files
//pcd_files.push_back(file2);
//pcd_files.push_back(file3);
pcd_files.push_back(file4);

vector<Plane> planes;
vector<Scene> scenes;

Orientation base_orientation;

Matrix3b3 R_del;

//Read in IE output and save to matrix
MatrixXd GNSS_INS_data;
Read_Mat("LiDAR_Georeferencing_Output_noheader.txt", GNSS_INS_data); //This is already done before the loop - is this necessary?
																	 //Read_Mat("Orientation.txt", GNSS_INS_data);
//Read in EOP's for each orientation
MatrixXd Orientation_EOP;
Read_Mat("C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Edmond_Cross_EOP2.txt", Orientation_EOP);
//Read_Mat("C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\AppleEOPs_noheader.txt", Orientation_EOP);


//for (int i = 0; i < pcd_files.size(); i++)
//{
//	clog << "\n-------------------------Starting on Scene " << i << "-------------------------------------------------------\n";
//
//	//Clear vector of planes
//	planes.clear();
//
//	Scene temp_scene;
//
//	// ---------------------------------------STEP 1: Load PCD Scene Data-----------------------------------------------------------------------------------------
//	//clog << "\n-------------------------STEP 1: Load PCD Scene Data-------------------------------------------------------\n";
//
//	//std::clog << "Opening file: " << pcd_files[i] << " (can take up to 5 minutes)" << endl;
//	clog << "Loading file....";
//	PointCloudXYZptr Novatel_cloud(new PointCloudXYZ);
//	if (!Read_Lidar_points(pcd_files[i], Novatel_cloud))
//	{
//		clog << "\n\nProgram will close\n";
//		return -1;
//	};
//
//
//	//clog << "\n-------------------------STEP 2: Filter Data-------------------------------------------------------\n";
//
//
//	// Create the filtering object and downsample. (USE SUBSAMPLING INSTEAD)
//	clog << "Filtering Dataset....";
//	filter_and_downsample(Novatel_cloud, 0.1f);
//
//
//	//clog << "\n-------------------------STEP 3: Fit all planes-----------------------------------------------------\n";
//
//	planes = FitPlanes(Novatel_cloud);
//
//
//
//	// Find the largest planes
//	std::sort(planes.begin(), planes.end(), sort_cloud); // sort based off cloud size
//
//	planes.resize(5); //truncate to keep largest planes
//
//
//
//	//clog << "\n-------------------------STEP 4: Downsample pts on Planes----------------------------------------------------\n";
//
//	clog << "\nDownsampling.....\n\n";
//
//	for (int i = 0; i < planes.size(); i++) {
//		remove_outliers(planes[i].points_on_plane, 100, 1);
//		filter_and_downsample(planes[i].points_on_plane, 0.5f);
//	}
//	//Remove Small Planes
//	int removed = 0;
//	for (int i = planes.size()-1; i >= 0; i--)
//	{
//		clog << "Cloud size: " << planes[i].points_on_plane->size() << endl;
//		if (planes[i].points_on_plane->size() < 250) {
//			planes.erase(planes.begin() + i); //Too few points
//			clog << "removed. Too small\n";
//		}
//	}
//
//		//save planes
//		save_planes(planes);
//
//
//		//clog << "\n-------------------------STEP 4: Get IE GNSS/INS OBS-----------------------------------------------------\n";
//
//
//		//GNSS
//		temp_scene.scene_orientation.X = Orientation_EOP(i,2);
//		temp_scene.scene_orientation.Y = Orientation_EOP(i, 3);
//		temp_scene.scene_orientation.Z = Orientation_EOP(i, 4);
//		//INS
//		temp_scene.scene_orientation.omega = Orientation_EOP(i, 5);
//		temp_scene.scene_orientation.phi = Orientation_EOP(i, 6);
//		temp_scene.scene_orientation.kappa = Orientation_EOP(i, 7);
//		
//
//		// DEBUG
//		cout << "EOP of ORIENTATION " << i << "\n\n";
//		cout << "\tX:\t" << temp_scene.scene_orientation.X << endl;
//		cout << "\tY:\t" << temp_scene.scene_orientation.Y << endl;
//		cout << "\tZ:\t" << temp_scene.scene_orientation.Z << endl;
//		cout << "\tOmega:\t" << temp_scene.scene_orientation.omega << endl;
//		cout << "\tPhi:\t" << temp_scene.scene_orientation.phi << endl;
//		cout << "\tKappa:\t" << temp_scene.scene_orientation.kappa << endl;
//
//
//		// VISUALIZE
//		//visualize_planes(planes);
//		
//		//Add to scene
//		temp_scene.planes = planes;
//		scenes.push_back(temp_scene);
//
//		//DEBUG
//		cout << "\n\nPlane equations for " << i << endl;
//		for (int k = 0; k < planes.size(); k++)
//		{
//			cout << "\t" << planes[k].a1 << "\t" << planes[k].a2 << "\t" << planes[k].a3 << "\t" << planes[k].b << endl;
//		}
//		cout << "-------------------------- END -------------------------\n\n";
//	}

	// DEBUG INPUT
	scenes = LoadDebugData();

	clog << "\n-------------------------Finished finding planes -------------------------------------------------------\n";
	clog << "Matching planes....";
	// TODO: MATCH PLANES
	UniquePlanes unique_planes = match_scenes(scenes);

	clog << "\n\nMapping vector PRE REMOVE:\n";
	print_vector(unique_planes.mapping_vec);
	clog << "Done.\nRemoving unfrequent planes....";
	//Remove less frequent planes. 
	int num_removed = remove_unfrequent(unique_planes);
	clog << "Done. Removed " << num_removed << " infrequent planes.\n";

	//DEBUG
	save_planes(unique_planes.unique_planes);



	// ------------DEBUG
	clog << "\nNumber of planes:\t" << unique_planes.unique_planes.size() << endl;
	clog << "\Frequency vector:\n";
	print_vector(unique_planes.frequency);
	clog << "\nMapping vector POST REMOVE:\n";
	print_vector(unique_planes.mapping_vec);

	// -----------------


	// Make output vectors
	vector<RowVectorXd> point_details; // X, Y, Z, 
	vector<RowVectorXd> scene_details; // X, Y, Z, Omega, Phi, Kappa
	vector<RowVectorXd> plane_details; // A1, A2, A3, d

	create_bundle_observations(scenes, unique_planes, point_details, scene_details, plane_details);

	char * scene_file = "Scene_Details.txt";
	char * plane_file = "Plane_Details.txt";
	char * point_file = "Point_Details.txt";

	clog << "\n---------------------------Scene---------------------\n";
	clog << "\tThere are " << scene_details.size() << " scenes\n";
	print_vector(scene_details);
	print_vector(scene_details, scene_file);

	clog << "\n---------------------------Plane---------------------\n";
	clog << "\tThere are " << plane_details.size() << " planes\n";
	print_vector(plane_details);
	print_vector(plane_details, plane_file);

	clog << "\n---------------------------Points---------------------\n";
	clog << "\tThere are " << point_details.size() << " points\n";
	print_vector(point_details);
	print_vector(point_details, point_file);

	
	
	// BUNDLE ADJUSTMENT
	MatrixXd planeparams;
	MatrixXd scanparams;
	MatrixXd pointparams;	
	
	//Convert vector of vectors to MatrixXd for plane, scene, and lidar point parameters
	vec2mat(plane_details, planeparams, 4);
	vec2mat(scene_details, scanparams, 6);
	vec2mat(point_details, pointparams, 3);
	
	//Get number of planes, scans, lidar points
	int numPlanes = planeparams.rows();
	int numScans = scanparams.rows();
	int numLidPts = pointparams.rows();
	
	//total number of unknowns
	int u = 6 + 4*numPlanes + 6*numScans; 	
	
	//Initial boresight parameters set to zero
	MatrixXd bs_params = MatrixXd::Zero(6,1);
		
	//Initialize matrices for LS
	MatrixXd A(numScans*6+numLidPts, u);
	MatrixXd A_full(numScans*6+numPlanes+numLidPts, u);
	MatrixXd N(0,0);	
	MatrixXd w(numScans*6+numLidPts, 1);
	MatrixXd w_full(numScans*6+numPlanes+numLidPts, 1); 
	MatrixXd U(0,0);	
	MatrixXd H(numPlanes, u);
	MatrixXd V(numPlanes ,1);	
	MatrixXd B(u+numPlanes, u+numPlanes);
	MatrixXd C(u+numPlanes,1);
	MatrixXd Y(u+numPlanes,1);
	
	//Set Clo using GNSS_INS_data - TO DO **************************************
	MatrixXd Clo = MatrixXd::Identity(numScans * 6 + numPlanes + numLidPts, numScans * 6 + numPlanes + numLidPts);
		
	//Compute P using Clo. Use apriorivf=1
	MatrixXd P = Clo.inverse();
	
	double mean_Y = 1;
	
	while(mean_Y > 0.000001)
	{
		computeAandw(A_full, A, H, w_full, w, V, u, numPlanes, numScans, numLidPts, bs_params, planeparams, scanparams, pointparams, GNSS_INS_data);
		//cout << "A=" << A_full << "\n\n";
		//cout << "w=" << w << "\n\n";

		N = A_full.transpose() * P * A_full; //Not sure if N is based on the full A, P with plane eqn derivatives or not???******************************
		/* cout << "A^T=" << A_full.transpose() << "\n\n";
		cout << "P" << P << "\n\n";
		cout << "N=" << N << "\n\n"; */

		U = -1 * A_full.transpose() * P * w_full; //Not sure if A, P, w here are based on the full A with plane eqn derivatives or not???******************************
		//cout << "U=" << U << "\n\n";

		B.block(0, 0, u, u) = N;
		B.block(u, 0, numPlanes, u) = H;
		B.block(0, u, u, numPlanes) = H.transpose();
		B.block(u, u, numPlanes, numPlanes) = MatrixXd::Constant(numPlanes, numPlanes, 0);
	
		C.block(0,0,u,1) = U;
		C.block(u,0,numPlanes,1) = V;
		
		Y = (B.transpose()*B).inverse() * B.transpose() * C;
		//cout << "Y=" << Y << "\n\n";

		//Calculate mean value of Y for termination condition
		mean_Y = 0;
		for(int i=0; i<(u+numPlanes); i++)
		{
			mean_Y = mean_Y + Y(i,0);
		}
		mean_Y = mean_Y/(u+numPlanes);
		//cout << "mean Y=" << mean_Y << "\n\n";

		//Update unknowns
		
		bs_params = bs_params + Y.block(0,0,6,1);

		int Yindex = 6;
		for (int i = 0; i < planeparams.rows(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				planeparams(i, j) = Y(Yindex,0);
				Yindex++;
			}
		}

		for (int i = 0; i < scanparams.rows(); i++)
		{
			for (int j = 0; j < 6; j++)
			{
				scanparams(i, j) = Y(Yindex, 0);
				Yindex++;
			}
		}

	}
	

	// GEOREFERENCE
	// TODO: Read in LiDAR frame -- or extract frame from full dataset? 
	MatrixXd lidar_data;
	//char * lidar_file = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\CalibrationAlgorithm\\BoresightAlgorithm\\Build\\Cross_01(Frame 0000).txt";
	//Read_Mat(lidar_file, lidar_data);

	// TODO: Match timestamps of GNSS+INS data to that of LiDAR data
	// Assumes that lidar_data is all points associated with one frame
	// Need to fix below function so that we can equate reference frame of timestamps
	//MatrixXd combined_data = merge_data(GNSS_INS_data, lidar_data);

	// TODO: Call georeferencing function
	//MatrixXd output_cloud = georeference_lidar_point(combined_data, boresight_leverarm, boresight_angles);

	// OPTIONAL: Convert final MatrixXd to pcd::PointCloud or a custom struct

	clog << "\n\n FINISHED CALIBRATION\n\n";


	return 0;
}