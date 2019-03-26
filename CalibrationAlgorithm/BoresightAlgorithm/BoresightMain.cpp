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


	fprintf(stdout, "    _________|T|_______         											\n");
	fprintf(stdout, "   |,-----.,-----.,---.\        											\n");
	fprintf(stdout, "   ||  IF || THE ||MAP \\       											\n");
	fprintf(stdout, "   |`-----'|-----||-----\`----. 											\n");
	fprintf(stdout, "   [       |    -||-   _|    (| 											\n");
	fprintf(stdout, "   [  ,--. |_FITS||___/.--.   | 											\n");
	fprintf(stdout, "   =-(( `))-----------(( `))-== 											\n");
	fprintf(stdout, "      `--'             `--'     											\n");
	
	
	fprintf(stdout, "\n---------------------------------------------------------------------------------------\n");
	
	//Initialize variables 
	vector<char *> pcd_files;
	//char * file1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation1_b.pcd";
	//char * file2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation2_b.pcd";
	//char * file3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation3_b.pcd";
	//char * file4 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Orientation4_b.pcd";
	
	/*char * file1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Orientation1_b.pcd";
	char * file2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Orientation2_b.pcd";
	char * file3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Orientation3_b.pcd";
	char * file4 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Orientation4_b.pcd";*/
	
	//char * file1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
	//char * file2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points_dec5.pcd";
	//char * file3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
	
	//char * file1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
	//char * file2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
	//char * file3 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\Data\\FirstDataset\\All_points.pcd";
	
	//pcd_files.push_back(file1); //vector of input files
	//pcd_files.push_back(file2);
	//pcd_files.push_back(file3);
	//pcd_files.push_back(file4);

	
	Orientation base_orientation;
	
	Matrix3b3 R_del;
	
	//Read in IE output and save to matrix
	MatrixXd GNSS_INS_data;
	Read_Mat("LiDAR_Georeferencing_Output_noheader.txt", GNSS_INS_data); //This is already done before the loop - is this necessary?
																		 //Read_Mat("Orientation.txt", GNSS_INS_data);
	//Read in EOP's for each orientation
	MatrixXd Orientation_EOP;
	//Read_Mat("C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Edmond_Cross_EOP2.txt", Orientation_EOP);
	//Read_Mat("C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\Apple_EOP_Edmond.txt", Orientation_EOP);

	Read_Mat("C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\CalibrationAlgorithm\\BoresightAlgorithm\\Build\\Edmond_Cross_EOP2.txt", Orientation_EOP);
	//Read_Mat("C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\AppleWarehouse\\AppleEOPs_noheader.txt", Orientation_EOP);

	// Load the files into scenes
	//vector<Scene> scenes = load_scenes(pcd_files, Orientation_EOP);

	//``````````````````````````````````````````````````````````````````````START COMMENTING: DEBUGGING GEOREFERENCING STUFF


	 //DEBUG INPUT
	vector<Scene> scenes = LoadDebugData();


	//Get the Shiftdown for distance measurement simplification
	vector<double> shiftdown;
	find_apply_shiftdown(scenes, shiftdown);
	


	clog << "\n-------------------------Finished finding planes -------------------------------------------------------\n";
	clog << "Matching planes....";
	// TODO: MATCH PLANES
	UniquePlanes unique_planes = match_scenes(scenes);
	clog << "Done.\nRemoving unfrequent planes....";

	clog << "\n\nMapping vector POST REMOVE:\n";
	print_vector(unique_planes.mapping_vec);
	//Remove less frequent planes. 
	int num_removed = remove_unfrequent(unique_planes);
	clog << "\n\nDone. Removed " << num_removed << " infrequent planes.\n";




	// ------------DEBUG
	clog << "\nNumber of planes:\t" << unique_planes.unique_planes.size() << endl;
	clog << "\Frequency vector:\n";
	print_vector(unique_planes.frequency);
	clog << "\nMapping vector POST REMOVE:\n";
	print_vector(unique_planes.mapping_vec);

	// -----------------


	// Make output vectors
	vector<RowVectorXd> point_details; // X, Y, Z, Unique Plane, Scene
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
	//print_vector(point_details);
	//print_vector(point_details, point_file);
	
	
	// BUNDLE ADJUSTMENT
	MatrixXd planeparams;
	MatrixXd scanparams;
	MatrixXd pointparams;
	MatrixXd GNSSINSparams = Orientation_EOP.block(0, 2, Orientation_EOP.rows(), 6);

	//Convert vector of vectors to MatrixXd for plane, scene, and lidar point parameters
	vec2mat(plane_details, planeparams, 4);
	vec2mat(scene_details, scanparams, 6);
	vec2mat(point_details, pointparams, 5);

	Write_Mat("C:\\Users\\kiera.fulton2\\Desktop\\esther\\planeparams.txt", planeparams, 5);
	Write_Mat("C:\\Users\\kiera.fulton2\\Desktop\\esther\\scanparams.txt", scanparams, 5);
	Write_Mat("C:\\Users\\kiera.fulton2\\Desktop\\esther\\pointparams.txt", pointparams, 5);

	//Get number of planes, scans, lidar points
	int numPlanes = planeparams.rows();
	int numScans = scanparams.rows();
	int numLidPts = pointparams.rows();

	//total number of unknowns
	int u = 6 + 4 * numPlanes + 6 * numScans;

	//Initial boresight parameters set to zero
	MatrixXd bs_params = MatrixXd::Zero(6, 1);

	//Initialize matrices for LS
	MatrixXd A(numScans * 6 + numLidPts, u);
	MatrixXd A_full(numScans * 6 + numPlanes + numLidPts, u);
	MatrixXd N(0, 0);
	MatrixXd w(numScans * 6 + numLidPts, 1);
	MatrixXd w_full(numScans * 6 + numPlanes + numLidPts, 1);
	MatrixXd U(0, 0);
	MatrixXd H(numPlanes, u);
	MatrixXd V(numPlanes, 1);
	MatrixXd B(u + numPlanes, u + numPlanes);
	MatrixXd C(u + numPlanes, 1);
	MatrixXd Y = MatrixXd::Ones(u + numPlanes, 1);

	//Change hardcoded values to actual values******************************************************
	//Set P. Clo is a diagonal matrix so inverse of Clo is a diagonal matrix with each diagonal element inversed.
	MatrixXd P_full = MatrixXd::Identity(numScans * 6 + numPlanes + numLidPts, numScans * 6 + numPlanes + numLidPts);
	for (int i = 0; i < (numScans * 6); i = i + 6)
	{
		P_full(i, i) = 1 / 0.003;
		P_full(i + 1, i + 1) = 1 / pow(0.003, 2);
		P_full(i + 2, i + 2) = 1 / pow(0.004, 2);
		P_full(i + 3, i + 3) = 1 / pow(0.00069 * 180 / PI, 2);
		P_full(i + 4, i + 4) = 1 / pow(0.00070 * 180 / PI, 2);
		P_full(i + 5, i + 5) = 1 / pow((0.00466 + 0.00289 + 0.00339) / 3 * 180 / PI, 2);
	}
	P_full.block(numScans * 6, numScans * 6, numPlanes, numPlanes) = MatrixXd::Identity(numPlanes, numPlanes) * 1 / sqrt(pow(0.003, 2) + pow(0.003, 2) + pow(0.004, 2));
	P_full.block(numScans * 6 + numPlanes, numScans * 6 + numPlanes, numLidPts, numLidPts) = MatrixXd::Identity(numLidPts, numLidPts) * 1 / 0.03;

	MatrixXd P(numScans * 6 + numLidPts, numScans * 6 + numLidPts);
	P.block(0, 0, numScans * 6, numScans * 6) = P_full.block(0, 0, numScans * 6, numScans * 6);
	P.block(numScans * 6, numScans * 6, numLidPts, numLidPts) = P_full.block(numScans * 6 + numPlanes, numScans * 6 + numPlanes, numLidPts, numLidPts);


	//cout << endl << P << endl << endl;


	int iter = 0;
	//double mean_Y = 1;
	//mean_Y > 0.000001
	//cout << "Y: " << Y.block(0, 0, 6, 1);

	//for (int iter = 0; iter < 10; iter++)
	while (Y(0, 0) > 0.0001 || Y(1, 0) > 0.0001 || Y(2, 0) > 0.0001 || Y(3, 0) > 0.0001 || Y(4, 0) > 0.0001 || Y(5, 0) > 0.0001)
	{
		iter++;

		computeAandw(A_full, A, H, w_full, w, V, u, numPlanes, numScans, numLidPts, bs_params, planeparams, scanparams, pointparams, GNSSINSparams);
		//cout << "A=" << A_full << "\n\n";
		cout << "w=" << w << "\n\n";

		//N = A_full.transpose() * P_full * A_full; //Not sure if N is based on the full A, P with plane eqn derivatives or not???******************************
		/*cout << "A^T=" << A_full.transpose() << "\n\n";
		cout << "P" << P << "\n\n";*/
		N = A.transpose() * P * A;
		//cout << endl << "N=" << N << "\n\n";

		//U = -1 * A_full.transpose() * P_full * w_full; //Not sure if A, P, w here are based on the full A with plane eqn derivatives or not???******************************
		//cout << "U=" << U << "\n\n";
		U = -1 * A.transpose() * P * w;

		B.block(0, 0, u, u) = N;
		B.block(u, 0, numPlanes, u) = H;
		B.block(0, u, u, numPlanes) = H.transpose();
		B.block(u, u, numPlanes, numPlanes) = MatrixXd::Zero(numPlanes, numPlanes);

		C.block(0, 0, u, 1) = U;
		C.block(u, 0, numPlanes, 1) = V;

		MatrixXd BTB = B.transpose()*B;
		FullPivLU<MatrixXd> B_decomp(BTB);
		cout << "Is BTB invertible? " << B_decomp.isInvertible() << endl;

		Y = (B.transpose()*B).inverse() * B.transpose() * C;
		//cout << "Y=" << Y << "\n\n";

		//Calculate mean value of Y for termination condition
		/*mean_Y = 0;
		for (int i = 0; i<(u + numPlanes); i++)
		{
		mean_Y = mean_Y + Y(i, 0);
		}
		mean_Y = mean_Y / (u + numPlanes);*/
		//cout << "mean Y=" << mean_Y << "\n\n";

		//Update unknowns

		bs_params = bs_params + Y.block(0, 0, 6, 1);

		int Yindex = 6;
		for (int i = 0; i < planeparams.rows(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				planeparams(i, j) = Y(Yindex, 0);
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

		cout << "Iteration: " << iter << endl;
		cout << "Delta: " << Y.block(0, 0, 6, 1) << endl;

	}

	//Calculate obs residuals, aposteriori variance factor and Cxhat
	MatrixXd v = A_full * Y.block(0, 0, u, 1) + w_full;
	MatrixXd apostvf = v.transpose() * P * v / (numScans * 6 + numPlanes + numLidPts - u);

	FullPivLU<MatrixXd> lu_decomp(N);
	cout << "Is N invertible? " << lu_decomp.isInvertible() << endl;
	//cout << N << endl;
	MatrixXd Cxhat = apostvf * N.inverse();

	//Output: boresight parameters bs_params and Cxhat
	cout << "Final boresight parameters:" << bs_params << endl;
	cout << Cxhat(0, 0) << endl;
	
	//TEMP FOR DEBUGGING

	//``````````````````````````````````````````````````````````````````````STOP COMMENTING: DEBUGGING GEOREFERENCING STUFF

	// GEOREFERENCE
	//Create output file
	ifstream outputFile;
	outputFile.open("Georeferenced Point Cloud.txt", ios::in);
	if (outputFile.fail())
	{
		cout << "Problem creating final output file. " << endl;
	}
	else
	{
		cout << "Created final output file..." << endl;
	}

	// TODO: Read in LiDAR frame -- or extract frame from full dataset? 
	MatrixXd lidar_data;
	MatrixXd current_lidar_point;
	current_lidar_point.resize(1,7); // or however many columns there are in the lidar data
	//Convert pcap of all point to text file
	char * lidar_file = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\IfTheMapFits\\CalibrationAlgorithm\\BoresightAlgorithm\\Build\\Cross_01(Frame 0000).txt";
	Read_Mat(lidar_file, lidar_data);

	//TODO: Find start time from IE file
	double start_time = GNSS_INS_data(0, 0);
	cout << "Start time = " << start_time << endl;

	//Convert LiDAR time to GPS week seconds
	double hour;
	int day;
	get_hour_day(start_time, &hour, &day);
	cout << "Hour: " << hour << " Day: " << day << endl;

	double prev_time = 0.0;

	for (int i = 1; i < size(lidar_data); i++){
		//Retrieve record in lidar data
		double current_lidar_time = lidar_data(i, 10);
		// Find time of lidar data point in GPS seconds

		if (prev_time > current_lidar_time)
		{
			hour++;
			if (hour > 23.0)
			{
				hour = 0;
				day++;
			}
		}

		double current_gpssec = (day * 86400) + (hour * 3600) + (current_lidar_time / (3.6*pow(10, 6)));
		
		
		//TODO: function to round time to the nearest quarter second
		current_gpssec = round_time(current_gpssec);
		cout << "current_gpssec (rounded): " << current_gpssec << endl;

		MatrixXd combined_data;
		combined_data = merge_data(GNSS_INS_data, lidar_data, current_gpssec);

		//TODO: Call georeferencing function
		MatrixXd boresight_leverarm, boresight_angles;
		boresight_leverarm << 0.0, 0.0, 0.0;
		boresight_angles << 180.0, 0.0, 90.0;
		MatrixXd output_cloud = georeference_lidar_point(combined_data, boresight_leverarm, boresight_angles);

		prev_time = current_lidar_time;

		//TODO: Read out single georeferenced point and append output file

	}


	clog << "\n\n FINISHED CALIBRATION\n\n";

	cin.get();
	return 0;
}