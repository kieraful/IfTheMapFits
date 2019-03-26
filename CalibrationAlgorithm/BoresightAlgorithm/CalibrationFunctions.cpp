#include "CalibrationFuntions.h"



bool Read_Lidar_points(char *filename, PointCloudXYZptr cloud)
{
	pcl::PCDReader reader;
	pcl::ScopeTime readfilescope("File Read");
	{
		if (reader.read(filename, *cloud) == -1) {
			std::cerr << "File could not be opened:\n\t" << filename << endl;
			return false;
		}
	}
	return true;
}

//Receives the name of a file "FileName" containing a numerical matrix, and read the matrix data into variable "m"
void Read_Mat(char *FileName, MatrixXd& m) {

	m.resize(0, 0);

	ifstream matfile;
	matfile.open(FileName, ios::in); //to open the file and start reading from the beginning
	if (matfile.fail()) //check if the file is opened successfully
	{
		cout << "\nThere was a problem reading the following file: " << endl << FileName << endl;
		//exit (EXIT_FAILURE);
		return;
	}
	else {
		cout << "\nFile read correctly. Continuing.\n";
	}

	char* readlinechr = new char[MaxMatSize];
	vector<double> v_all;
	int nrow = 0;

	while (matfile.getline(readlinechr, MaxMatSize, '\n')) {
		nrow++;
		int stln = strlen(readlinechr);
		char* readlinestr = new char[stln + 1];
		for (int i = 0; i<stln; i++)
		{
			readlinestr[i] = readlinechr[i];
		}

		readlinestr[stln] = '\0';

		stringstream rowstream(readlinestr);
		double value;
		while (!rowstream.eof()) {
			rowstream >> value;
			v_all.push_back(value);
		}
	}
	matfile.close();

	int ncol = v_all.size() / nrow;
	m.resize(nrow, ncol);

	for (int i = 0; i<nrow; i++) {
		for (int j = 0; j<ncol; j++) {
			m(i, j) = v_all.at(i*ncol + j);
		}
	}

	return;
};

//Recives the name of a file "FileName", and writes the data in the matrix "m" to this file, with fixed precision "decimal_precision"
void Write_Mat(char *FileName, MatrixXd &m, int decimal_precision) {

	ofstream matfile;
	matfile.open(FileName, ios::out); //to open the file and start writing from the beginning (Notice! over writing)
	if (matfile.fail()) //check if the file is opened successfully
	{
		cout << "There was a problem reading the follwoing file: " << endl << FileName << endl;
		//exit (EXIT_FAILURE);
		return;
	}
	matfile.flags(ios::fixed);
	matfile.precision(decimal_precision);
	matfile << m;
	matfile.close();
	return;
};

//Receives three rotation angles (Omega,Phi,Kappa) and returns the rotation matrix from the object to the image space (Rot_g2i)
//Note that the rotation from image to object space is the transpose of "Rot_g2i"
void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i) {
	Matrix3b3 Mw0;
	Matrix3b3 Mf0;
	Matrix3b3 Mk0;

	// compute R_g_to_i and return it to the Rot_g2i

	Mw0 << 1, 0, 0,
		0, cosd(Omega), sind(Omega),
		0, -sind(Omega), cosd(Omega);

	Mf0 << cosd(Phi), 0, -sind(Phi),
		0, 1, 0,
		sind(Phi), 0, cosd(Phi);

	Mk0 << cosd(Kappa), sind(Kappa), 0,
		-sind(Kappa), cosd(Kappa), 0,
		0, 0, 1;

	Rot_g2i = Mk0 * Mf0*Mw0;
};

void Convert_R_to_Angles(Matrix3b3 R, double& Omega, double& Phi, double& Kappa) {

	double A11 = R(0, 0);
	double A12 = R(0, 1);
	double A13 = R(0, 2);
	double A21 = R(1, 0);
	double A22 = R(1, 1);
	double A23 = R(1, 2);
	double A33 = R(2, 2);

	if (A13 != 1 && A13 != -1) {

		Phi = asin(A13);
		Omega = atan2(-A23 / cos(Phi), A33 / cos(Phi));
		Kappa = atan2(-A12 / cos(Phi), A11 / cos(Phi));

	}

	else if (A13 == 1) {
		Phi = PI / 2;
		Omega = 0; //arbitrary
		Kappa = -Omega + atan2(A21, A22);

	}

	else if (A13 == -1) {
		Phi = -PI / 2;
		Omega = 0;
		Kappa = Omega + atan2(A21, A22);
	}

	return;
};

//Given corresponding points:
//       xy_i1: x in pixels, y in pixels (image observations of points on image 1)
//		 xy_i2: x in pixels, y in pixels (image observations of the same points on image 2)
// So, the size of "xy_i1" and "xy_i2" must be the same
//Finds the conditioning homographies, H1 and H2
void  Normalization_Condition(MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& H1, MatrixXd& H2) {

	int n_p = xy_i1.rows();

	double cx1 = (xy_i1.block(0, 0, n_p, 1)).array().mean();
	double cy1 = (xy_i1.block(0, 1, n_p, 1)).array().mean();
	double cx2 = (xy_i2.block(0, 0, n_p, 1)).array().mean();
	double cy2 = (xy_i2.block(0, 1, n_p, 1)).array().mean();

	MatrixXd dx1 = xy_i1.block(0, 0, n_p, 1) - cx1 * (MatrixXd::Ones(n_p, 1));
	MatrixXd dy1 = xy_i1.block(0, 1, n_p, 1) - cy1 * (MatrixXd::Ones(n_p, 1));
	MatrixXd dx2 = xy_i2.block(0, 0, n_p, 1) - cx2 * (MatrixXd::Ones(n_p, 1));
	MatrixXd dy2 = xy_i2.block(0, 1, n_p, 1) - cy2 * (MatrixXd::Ones(n_p, 1));

	MatrixXd Temp = (dx1.array().square() + dy1.array().square()).array().sqrt();
	double d1 = Temp.array().mean();

	Temp = (dx2.array().square() + dy2.array().square()).array().sqrt();
	double d2 = Temp.array().mean();

	H1.setZero(3, 3);
	H2.setZero(3, 3);

	H1 << sqrt(2.0) / d1, 0, -(sqrt(2.0) / d1 * cx1),
		0, sqrt(2.0) / d1, -(sqrt(2.0) / d1 * cy1),
		0, 0, 1;
	H2 << sqrt(2.0) / d2, 0, -(sqrt(2.0) / d2 * cx2),
		0, sqrt(2.0) / d2, -(sqrt(2.0) / d2 * cy2),
		0, 0, 1;


	return;
}
void Find_closest_points(int num_find_points, LidarPt point, vector<LidarPt> search_points, vector<LidarPt>& closest_pts)
{
	/*
		Code to fine the closest points to a given point. 
		inputs:
			num_find_points: (`int`) The number of points to find that are closest
			point: (`struct`) The start point
			search_points: (`vector<lidar>`) The search space
	
	
	*/
	LidarPt best_point;
	for (int j = 0; j < num_find_points; j++) {

		int num_search_points = search_points.size();
		double best_index = -1;
		double best_dist = -1;
		double dist;

		for (int i = 0; i < num_search_points; i++) {

			// LINEAR SEARCH
			dist = euclidian_dist(point.x, point.y, point.z, search_points[i].x, search_points[i].y, search_points[i].z);
			if (dist < best_dist) {
				best_dist = dist;
				best_index = i;

			}

		}

		best_point.x = search_points[best_index].x;
		best_point.y = search_points[best_index].y;
		best_point.z = search_points[best_index].z;

		//Remove the row of the best point from search points
		search_points.erase(search_points.begin() + best_index);
		// Add point to closest points
		closest_pts.push_back(best_point);

	}



}
double euclidian_dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}




vector<Plane> FitPlanes(PointCloudXYZptr in_cloud, int max_planes, bool make_files) {

	/*
	pcl::PointCloud<pcl::PointXYZ> cloud_filtered:		 This is the filtered point cloud in which planes of interest lie
	int max_planes:										 This is the maximum number of planes to find in the cloud_filtered
	
	*/

	double percent_cloud = 30;

	if (max_planes > 0)
	{
		clog << " Fitting planes. Will stop at " << percent_cloud << "% remaining cloud or " << max_planes << " planes.\n\n";
	}
	else
	{
		clog << " Fitting planes. Will stop at " << percent_cloud << "% remaining cloud\n\n";
		max_planes = std::numeric_limits<int>::max();
	}
	//Initializers
	int n_planes = 0; // Number of planes found in dataset, init to 0
	vector<Plane> planes;
	Plane temp_plane;
	pcl::PCDWriter writer; //writer object for point clouds
	double plane_buffer = 25.0*(PI / 180.0f);//Degree offset from search plane to allow
	double size, remaining_p, remaining_pts;

	// Fill in the cloud data
	pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
	

	// Create the segmentation object
	pcl::SACSegmentation<pcl::PointXYZ> seg;

	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PARALLEL_PLANE); //only want points perpendicular to a given axis
	seg.setMaxIterations(1000);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setDistanceThreshold(0.15); // keep points within 0.15 m of the plane
	Eigen::Vector3f axis = Eigen::Vector3f(0.0, 0.0, 1.0); //z axis
	seg.setAxis(axis);
	seg.setEpsAngle(plane_buffer); // plane can be within 30 degrees of X-Z plane
	

	// Create the filtering object
	pcl::ExtractIndices<pcl::PointXYZ> extracter;

	int i = 0, nr_points = (int)in_cloud->points.size();
	// While 30% of the original cloud is still there
	while (in_cloud->points.size() > (percent_cloud/100) * nr_points && n_planes < max_planes)
	{
		// Initialize clouds
		PointCloudXYZptr cloud_p(new PointCloudXYZ), cloud_f(new PointCloudXYZ), cloud_temp;
		//Alert user to plane fitting

		// Segment the largest planar component from the remaining cloud
		seg.setInputCloud(in_cloud);
		seg.segment(*inliers, *coefficients);
		if (inliers->indices.size() == 0)
		{
			cerr << "\n\n\tERROR: No planes found in dataset!\n\n" << endl;
			break;
		}

		n_planes++;

		remaining_p = (1 - (double)in_cloud->points.size() / (double)nr_points)*100; //remaining percent
		remaining_pts = nr_points - in_cloud->points.size(); //remaining points

		fprintf(stdout, "\tSuccessfully fitted %0.3f %% of cloud. %0.1f points and %i planes\r", remaining_p, remaining_pts, n_planes);


		//Find plane parameters
		temp_plane.a1 = coefficients->values[0];
		temp_plane.a2 = coefficients->values[1];
		temp_plane.a3 = coefficients->values[2];
		temp_plane.b = coefficients->values[3];

		// Extract the plane inliers
		extracter.setInputCloud(in_cloud);
		extracter.setIndices(inliers);
		extracter.setNegative(false);
		extracter.filter(*cloud_p);

		//Save to plane struct
		temp_plane.points_on_plane = cloud_p;

		//get size
		size = temp_plane.points_on_plane->size();

		//pusback vector of planes
		planes.push_back(temp_plane);


		//write plane to file
		if (make_files)
		{
			std::stringstream ss;
			ss << "Cloud_Plane_" << i << ".pcd";
			writer.write<pcl::PointXYZ>(ss.str(), *cloud_p, false);
		}
		

		// Create the filtering object
		extracter.setNegative(true);
		extracter.filter(*cloud_f);
		in_cloud.swap(cloud_f);
		i++;

	}

	return planes;

}

void remove_outliers(PointCloudXYZptr &input_cloud, double search_n, double std_mult)
{

	// Create the filtering object
	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> out_remove;
	out_remove.setInputCloud(input_cloud);
	out_remove.setMeanK(search_n);
	out_remove.setStddevMulThresh(std_mult);
	out_remove.filter(*input_cloud);

}

void filter_and_downsample(PointCloudXYZptr &input_cloud, float leaf_size)
{


	PointCloudXYZptr filtered_cloud(new PointCloudXYZ);
	pcl::VoxelGrid<pcl::PointXYZ> vox_grid;
	pcl::ScopeTime filterscope("Filtering dataset");
	{
		vox_grid.setInputCloud(input_cloud);
		vox_grid.setLeafSize(leaf_size, leaf_size, leaf_size);
		vox_grid.filter(*input_cloud);
	}

}

bool sort_cloud(Plane plane_1, Plane plane_2)
{
	int num_1 = plane_1.points_on_plane->size();
	int num_2 = plane_2.points_on_plane->size();
	return (num_1 > num_2);
}

bool sort_planes(Vector3d vec_1, Vector3d vec_2)
{
	int num_1 = vec_1[0];
	int num_2 = vec_2[0];
	return (num_1 > num_2);
}

void visualize_planes(vector<Plane> planes)
{
	clog << "Visualizing " << planes.size() << " clouds with planes...\n";
	pcl::visualization::PCLVisualizer viewer("If The Map Fits");
	pcl::ModelCoefficients::Ptr plane(new pcl::ModelCoefficients);

	double r, g, b;
	std::stringstream ss;

	//For each plane 
	for (int i = 0; i < planes.size(); i++)
	{
		plane->values.resize(4);
		plane->values[0] = planes[i].a1;
		plane->values[1] = planes[i].a2;
		plane->values[2] = planes[i].a3;
		plane->values[3] = planes[i].b;


		//Add the plane
		ss.clear();
		ss << "Cloud_Plane_" << i;
		viewer.addPlane(*plane, ss.str(), 0);

		//Set the plane visuals

		r = ((double)rand() / (RAND_MAX));
		g = ((double)rand() / (RAND_MAX));
		b = ((double)rand() / (RAND_MAX));

		viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b /*R,G,B*/, ss.str(), 0);
		viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.6, ss.str(), 0);
		viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME, ss.str(), 0);

		//Add the point cloud
		ss.clear();
		ss << "Cloud_" << i;
		viewer.addPointCloud<pcl::PointXYZ>(planes[i].points_on_plane, ss.str());
	}


	viewer.addCoordinateSystem(0.5, "axis", 0);
	viewer.setBackgroundColor(0.05, 0.05, 0.05, 0);
	viewer.setPosition(800, 400);


	while (!viewer.wasStopped())
	{
		viewer.spinOnce();
	}
}


void visualize_cloud(PointCloudXYZptr cloud)
{

	clog << "Visualizing cloud...\n";
	pcl::visualization::PCLVisualizer viewer("If The Map Fits");

	while (!viewer.wasStopped())
	{
		viewer.addPointCloud<pcl::PointXYZ>(cloud);
		viewer.spinOnce();
	}



}

void save_plane(Plane save_plane, int identifyer)
{
	pcl::PCDWriter writer;
	std::stringstream ss;
	ss << "Cloud_Plane_" << identifyer << ".pcd";
	writer.write<pcl::PointXYZ>(ss.str(), *save_plane.points_on_plane, false);
}

void save_planes(vector<Plane> planes)
{

	//pcl::PCDWriter writer; //writer object for point clouds

	for (int i = 0; i < planes.size(); i++)
	{
		//std::stringstream ss;
		//ss << "Cloud_Plane_" << i << ".pcd";
		//writer.write<pcl::PointXYZ>(ss.str(), *planes[i].points_on_plane, false);
		save_plane(planes[i], i);
	}

}

MatrixXd georeference_lidar_point(MatrixXd data, MatrixXd boresight_LA, MatrixXd boresight_angles)
{
	//Vectors and matrices in georeferencing equation
	//1_2_3:
	//1: r or R (vector or rotation matrix)
	//2: subscript
	//3: superscript
	MatrixXd r_b_geo; //GNSS measurements -> geo coordinates of IMU center
	Matrix3b3 R_b_geo; //IMU measurements -> rotation of IMU in geo frame
	MatrixXd r_lidar_b; //Calibration lever arm
	Matrix3b3 R_lidar_b; //Calibration angles
	MatrixXd r_p_lidar; //LiDAR measurements
	MatrixXd r_p_geo; //ADDED TO OUTPUT -- Geo coords of the lidar point
	MatrixXd output; //Output: [timestamp, X_lidar, Y_lidar, Z_lidar]

	int num_points = data.rows();

	//Define dimensions
	r_b_geo.resize(3, 1); 
	R_b_geo.resize(3, 3);
	r_lidar_b.resize(3, 1);
	R_lidar_b.resize(3, 3);
	r_p_lidar.resize(3, 1);
	r_p_geo.resize(3, 1);
	output.resize(num_points, 5);


	double timestamp = 0.0;
	double vert_angle = 0.0;
	double horiz_angle = 0.0;
	double range = 0.0;

	//Iterating for each epoch
	for (int i = 0; i < num_points; i++)
	{
		timestamp = data(i, 0);

		//Populate matrices with applicable values
		r_b_geo << data(i, 5),
			data(i, 6),
			data(i, 7);

		Rotation_g2i(data(i, 8), data(i, 9), data(i, 10), R_b_geo); // Inputs: roll, pitch, azimuth, output matrix

		r_lidar_b << boresight_LA(0, 0),
			boresight_LA(0, 1),
			boresight_LA(0, 2); // Assumes boresight_LA read in is a 1x3 matrix 

		//R_lidar_b << boresight_angles(0, 0),
		//	boresight_angles(0, 1),
		//	boresight_angles(0, 2); // Assumes boresight_angles is a 1x3 matrix -> read these into Rotation_g2i matrix

		Rotation_g2i(boresight_angles(0, 0), boresight_angles(0, 1), boresight_angles(0, 2), R_lidar_b); // Assumes boresight_angles is a 1x3 matrix
			
		vert_angle = data(i, 2);
		horiz_angle = data(i, 3);
		range = data(i, 4);

		r_p_lidar << range * cos(vert_angle) * sin(horiz_angle),
			range * cos(vert_angle) * cos(horiz_angle),
			range*sin(vert_angle);

		r_p_geo = r_b_geo + (R_b_geo * r_lidar_b) + (R_b_geo * R_lidar_b * r_p_lidar);

		// push r_p_geo on to output with timestamp
		output(i, 0) = timestamp;
		output(i, 1) = data(i, 0); //Point ID
		output(i, 2) = r_p_geo(0, 0); //X_geo
		output(i, 3) = r_p_geo(1, 0); //Y_geo
		output(i, 4) = r_p_geo(2, 0); //Z_geo

	}

	// Output matrix will be: [timestamp, X, Y, Z]
	return output;
}

UniquePlanes match_scenes(vector<Scene> scenes)
{

	// For each scene, match planes. First scene is taken as base. 
	Scene base_scene;
	UniquePlanes unique;
	double del_omega, del_phi, del_kappa; // Defined as scene(i) - base
	Matrix3b3 R_del;
	RowVector3d target_plane_vec, base_plane_vec, target_rot_vec, mapping_temp, global_translation;
	vector<int> candidates;
	double thresh_orientation = 0.1;
	double best_dist, dot_prod, dist_temp, plane_dist, best_prod, az_test;
	int best_plane;
	Plane temp_plane;
	Orientation temp_orientation;

	//Add base scene to unique planes, as there is nothing to match yet
	for (int i = 0; i < scenes[0].planes.size(); i++) {


		plane_to_global(scenes[0].planes[i], scenes[0].scene_orientation);

		unique.unique_planes.push_back(scenes[0].planes[i]);
		mapping_temp << i , 0 , i;

		unique.mapping_vec.push_back(mapping_temp);
		unique.reference_orientations.push_back(scenes[0].scene_orientation);
		unique.frequency.push_back(1);

	}

	for (int i = 1; i < scenes.size(); i++) // For each Target scene i
	{
		for (int j = 0; j < scenes[i].planes.size(); j++) // for each plane j in target scene i
		{
			// Change plane orientation to global frame. 
			plane_to_global(scenes[i].planes[j], scenes[i].scene_orientation);
			// Target plane to vector
			target_plane_vec << scenes[i].planes[j].a1, scenes[i].planes[j].a2, scenes[i].planes[j].a3;
			
			candidates.clear();
			for (int k = 0; k < unique.unique_planes.size();k++) // for each unique plane k
			{
				
				// Base plane to matrix
				base_plane_vec << unique.unique_planes[k].a1, unique.unique_planes[k].a2, unique.unique_planes[k].a3;
				// ------- Dot Product Check -------------
				dot_prod = base_plane_vec.dot(target_plane_vec);
				if ((1 - abs(dot_prod)) < thresh_orientation)
				{
					candidates.push_back(k); // place candidate unique plane index in vector
				}

			}

			// If there is more than 1 candidate, check the relative distances

			// reset distance threshold
			best_dist = 10; // Planes should definitely not be more than 40 meters away from each other
			best_plane = -1;
			// For each candidate plane, find closest matching plane in base (Euclidian distance)
			if (candidates.size() > 1)
			{
				for (int m = 0; m < candidates.size(); m++)
				{
					//dist_temp = check_plane_dists(unique.reference_orientations[candidates[m]], scenes[i].scene_orientation, unique.unique_planes[candidates[m]], scenes[i].planes[j]);
					dist_temp = abs(abs(unique.unique_planes[candidates[m]].b) - abs(scenes[i].planes[j].b));


					//Azimuth check
					az_test = check_plane_az(unique.reference_orientations[candidates[m]], scenes[i].scene_orientation, unique.unique_planes[candidates[m]], scenes[i].planes[j]);



					//dist_temp = max(abs(unique.unique_planes[candidates[m]].b), abs(scenes[i].planes[j].b)) + 

					if (dist_temp < best_dist && az_test < 50)
					{
						//This is the best plane so far
						best_dist = abs(dist_temp);
						best_plane = candidates[m]; // save the base plane index
					}
				}
			}

			else if (candidates.size() == 1)
			{
				best_plane = candidates[0];
			}

			if (best_plane != -1)
			{
				// Unique plane match. This plane already exists. Add frequency
				// Add the best plane to the mapping matrix
				mapping_temp << best_plane, i, j;
				unique.mapping_vec.push_back(mapping_temp);
				//update_frequency
				unique.frequency[best_plane] = unique.frequency[best_plane] + 1;
				
			}
			else
			{
				// This plane is new. Add it to the unique planes.
				unique.unique_planes.push_back(scenes[i].planes[j]);
				mapping_temp << unique.unique_planes.size() - 1, i, j; // will be referencing the next unique plane
				unique.mapping_vec.push_back(mapping_temp);
				unique.frequency.push_back(1);
				unique.reference_orientations.push_back(scenes[i].scene_orientation);
			}

		}
		
	}

	return unique;
}

int remove_unfrequent(UniquePlanes & unique, int threshold)
{
	// Removes the less frequent planes from the mapping vector
	int removed = 0;

	for (int i = 0; i < unique.frequency.size(); i++)
	{
		if (unique.frequency[i] < threshold)
		{
			//This plane has too few referencing scenes. Remove from mapping vector
			removed++;
			for (int j = 0; j < unique.mapping_vec.size(); j++)
			{
				if ((int)unique.mapping_vec[j](0) == i)
				{
					unique.mapping_vec.erase(unique.mapping_vec.begin() + j);
				}
			}
		}
	}

	return removed;
}

void print_vector(vector<RowVector3d> print_vector)
{
	for (int i = 0; i < print_vector.size(); i++)
	{

		cout << "\t" << print_vector[i](0) << "\t" << print_vector[i](1) << "\t" << print_vector[i](2) << endl;

	}

}

void print_vector(vector<RowVectorXd> print_vector, char *filename)
{
	FILE *fp;
	fp = fopen(filename, "w");

	for (int i = 0; i < print_vector.size(); i++)
	{
		for (int j = 0; j < print_vector[i].cols(); j++)
		{
			fprintf(fp, "\t %0.3f ", print_vector[i](j));
		}
		fprintf(fp, "\n");

	}
}

void print_matrix(MatrixXd print_mat)
{

	for (int i = 0; i < print_mat.rows(); i++)
	{
		for (int j = 0; j < print_mat.cols(); j++)
		{
			printf("\t %0.3f ", print_mat(i,j));
		}
		printf("\n");

	}
}

void print_vector(vector<RowVectorXd> print_vector)
{
	for (int i = 0; i < print_vector.size(); i++)
	{
		for (int j = 0; j < print_vector[i].cols(); j++)
		{
			printf("\t %0.3f ",print_vector[i](j));
		}
		cout << endl;

	}
}

void print_vector(vector<int> print_vector)
{
	for (int i = 0; i < print_vector.size(); i++)
	{

		cout << "\t" << print_vector[i] << endl;

	}

}


double check_plane_dists(Orientation orient_base, Orientation orient_target, Plane p_base, Plane p_target)
{
	
	//Change to global distance. 
	Vector3d shiftdown, shifted_base, shifted_target, target_rot_vec, base_rot_vec;
	double d1, d2;
	//Find the shiftdown between EOP
	shiftdown(0) = orient_base.X - (orient_base.X - orient_target.X) + 1;
	shiftdown(1) = orient_base.Y - (orient_base.Y - orient_target.Y) + 1;
	shiftdown(2) = orient_base.Z - (orient_base.Z - orient_target.Z) + 1;

	shifted_base << orient_base.X - shiftdown(0), orient_base.Y - shiftdown(1), orient_base.Z - shiftdown(2);
	shifted_target << orient_target.X - shiftdown(0), orient_target.Y - shiftdown(1), orient_target.Z - shiftdown(2);


	//Find the plane equation vectors
	base_rot_vec << p_base.a1, p_base.a2, p_base.a3;
	target_rot_vec << p_target.a1, p_target.a2, p_target.a3;

	d1 = -1 * base_rot_vec.transpose() * shifted_base + p_base.b;
	d2 = -1 * target_rot_vec.transpose() * shifted_target + p_target.b;


	return abs(abs(d2) - abs(d1));

}

void find_apply_shiftdown(vector<Scene> &scenes, vector<double> &shiftdown)
{
	double average_X = 0;
	double average_Y = 0;
	double average_Z = 0;
	for (int i = 0; i < scenes.size(); i++)
	{
		average_X = average_X + scenes[i].scene_orientation.X;
		average_Y = average_Y + scenes[i].scene_orientation.Y;
		average_Z = average_Z + scenes[i].scene_orientation.Z;
	}

	average_X = (average_X) / scenes.size();
	average_Y = (average_Y) / scenes.size();
	average_Z = (average_Z) / scenes.size();

	// Find shiftdown values
	shiftdown.push_back(average_X);
	shiftdown.push_back(average_Y);
	shiftdown.push_back(average_Z);

	// Apply shiftdown values
	for (int j = 0; j < scenes.size(); j++)
	{
		scenes[j].scene_orientation.X = scenes[j].scene_orientation.X - shiftdown[0];
		scenes[j].scene_orientation.Y = scenes[j].scene_orientation.Y - shiftdown[1];
		scenes[j].scene_orientation.Z = scenes[j].scene_orientation.Z - shiftdown[2];
	}

	
}

//double check_plane_dists(Orientation orient_base, Orientation orient_target, Plane plane_base, Plane plane_target)
//{
//	// Checks the distance between two planes to see if they are in fact the same plane
//	// Done by finding circle intersection, and if it satisfies the base plane equation. 
//
//	bool result = false;
//	double d_base, d_target;
//	double x1, x2, y1, y2,z1, z2, r1, r2, d, a, h, px, py, int1x, int1y, int2x, int2y, check1, check2;
//
//
//	d_base = abs(plane_base.b); // Distance from origin in base scene
//	d_target = abs(plane_target.b); // Distance from origin in target scene
//
//	//find intersection points
//	x1 = orient_base.X;
//	y1 = orient_base.Y;
//	x2 = orient_target.X;
//	y2 = orient_target.X;
//	z1 = orient_base.Z;
//	z2 = orient_target.Z;
//
//	r1 = abs(plane_base.b);
//	r2 = abs(plane_target.b);
//
//	d = sqrt(pow((x1 - x2), 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
//
//	if (d > (r1 + r2))
//	{
//		//Does not intersect. Not the right plane. 
//		check1 = 100;
//		check2 = 100;
//
//	}
//	else if (d + min(r1, r2) <= max(r1, r2))
//	{
//
//		// One is in other
//		//TODO: Find closest point
//
//		check1 = abs(r1 - r2 - d);
//		check2 = 100;
//
//
//
//	}
//
//	else
//	{
//		a = (pow(r1, 2) - pow(r2, 2) + pow(d, 2)) / (2 * d);
//		h = sqrt(pow(r1, 2) - pow(a, 2));
//
//		px = x1 + (a*(x2 - x1)) / d;
//		py = y1 + (a*(y2 - y1)) / d;
//
//		int1x = px + (h*(y2 - y1)) / d;
//		int1y = py - (h*(x2 - x1)) / d;
//
//		int2x = px - (h*(y2 - y1)) / d;
//		int2y = py + (h*(x2 - x1)) / d;
//
//
//		// Now we have the intersection. Check if either one lies on the plane. 
//		check1 = abs(plane_base.a1*int1x + plane_base.a2*int1y + plane_base.a3*orient_target.Z - plane_base.b);
//		check2 = abs(plane_base.a1*int2x + plane_base.a2*int2x + plane_base.a3*orient_target.Z - plane_base.b);
//
//	}
//
//	return min(check1, check2);
//}

void rotate_scene(Scene & scene_target, Matrix3b3 R)
{
	//Rotate entire scene
	// Rotate all planes in scene
	// Rotate all Points in each plane
	// change Orientation to match base scene. 
	// This function does not move the location of the scene, obviously. Just the three orientation angles Omega, Phi, and Kappa
	// These values are taken from the rotation matrix R passed in. 

	Plane temp_plane;
	RowVector3d plane_vec;

	pcl::ScopeTime filterscope("Rotating scene");
	{
		//Rotate all planes in the scene. 
		for (int i = 0; i < scene_target.planes.size(); i++) // Each plane
		{
			plane_vec << scene_target.planes[i].a1, scene_target.planes[i].a2, scene_target.planes[i].a3;
			plane_vec = plane_vec * R;
			scene_target.planes[i].a1 = plane_vec(0);
			scene_target.planes[i].a2 = plane_vec(1);
			scene_target.planes[i].a3 = plane_vec(2);

			//Rotate all point in the plane
			for (int j = 0; j < scene_target.planes[i].points_on_plane->points.size(); j++) //each point on that plane
			{
				plane_vec << scene_target.planes[i].points_on_plane->points[i].x, scene_target.planes[i].points_on_plane->points[i].y, scene_target.planes[i].points_on_plane->points[i].z;
				plane_vec = plane_vec * R;
				scene_target.planes[i].points_on_plane->points[i].x = plane_vec(0);
				scene_target.planes[i].points_on_plane->points[i].y = plane_vec(1);
				scene_target.planes[i].points_on_plane->points[i].z = plane_vec(2);
			}

		}

	}
	//Get rotation angles. Change the scene orientation angles (Omega, Phi, Kappa)
	Convert_R_to_Angles(R, scene_target.scene_orientation.omega, scene_target.scene_orientation.phi, scene_target.scene_orientation.kappa);

}

void rotate_plane_points(Plane & plane_target, Matrix3b3 R)
{
	Vector3d plane_vec;
	//Rotate all point in the plane
	for (int j = 0; j < plane_target.points_on_plane->points.size(); j++) //each point on that plane
	{
		plane_vec << plane_target.points_on_plane->points[j].x, plane_target.points_on_plane->points[j].y, plane_target.points_on_plane->points[j].z;
		plane_vec = plane_vec.transpose() * R;
		plane_target.points_on_plane->points[j].x = plane_vec(0);
		plane_target.points_on_plane->points[j].y = plane_vec(1);
		plane_target.points_on_plane->points[j].z = plane_vec(2);
	}

}

vector<Scene> load_scenes(vector<char*> pcd_files, MatrixXd Orientation_EOP)
{

	vector<Plane> planes;
	vector<Scene> scenes;
	for (int i = 0; i < pcd_files.size(); i++)
	{
		clog << "\n-------------------------Starting on Scene " << i << "-------------------------------------------------------\n";

		//Clear vector of planes
		planes.clear();

		Scene temp_scene;

		// ---------------------------------------STEP 1: Load PCD Scene Data-----------------------------------------------------------------------------------------
		//clog << "\n-------------------------STEP 1: Load PCD Scene Data-------------------------------------------------------\n";

		//std::clog << "Opening file: " << pcd_files[i] << " (can take up to 5 minutes)" << endl;
		clog << "Loading file....";
		PointCloudXYZptr Novatel_cloud(new PointCloudXYZ);
		if (!Read_Lidar_points(pcd_files[i], Novatel_cloud))
		{
			clog << "\n\nCheck file " << i << "!!!\n";
		};


		//clog << "\n-------------------------STEP 2: Filter Data-------------------------------------------------------\n";


		// Create the filtering object and downsample. (USE SUBSAMPLING INSTEAD)
		clog << "Filtering Dataset....";
		filter_and_downsample(Novatel_cloud, 0.1f);


		//clog << "\n-------------------------STEP 3: Fit all planes-----------------------------------------------------\n";

		planes = FitPlanes(Novatel_cloud);



		// Find the largest planes
		std::sort(planes.begin(), planes.end(), sort_cloud); // sort based off cloud size

		planes.resize(5); //truncate to keep largest planes



						  //clog << "\n-------------------------STEP 4: Downsample pts on Planes----------------------------------------------------\n";

		clog << "\nDownsampling.....\n\n";

		for (int i = 0; i < planes.size(); i++)
		{
			remove_outliers(planes[i].points_on_plane, 100, 1);
			filter_and_downsample(planes[i].points_on_plane, 0.5f);
		}
		//Remove Small Planes
		int removed = 0;
		for (int i = planes.size() - 1; i >= 0; i--)
		{
			clog << "Cloud size: " << planes[i].points_on_plane->size() << endl;
			if (planes[i].points_on_plane->size() < 250)
			{
				planes.erase(planes.begin() + i); //Too few points
				clog << "removed. Too small\n";
			}
		}

		//save planes
		save_planes(planes);


		//clog << "\n-------------------------STEP 4: Get IE GNSS/INS OBS-----------------------------------------------------\n";


		//GNSS
		temp_scene.scene_orientation.X = Orientation_EOP(i, 2);
		temp_scene.scene_orientation.Y = Orientation_EOP(i, 3);
		temp_scene.scene_orientation.Z = Orientation_EOP(i, 4);
		//INS
		temp_scene.scene_orientation.omega = Orientation_EOP(i, 5);
		temp_scene.scene_orientation.phi = Orientation_EOP(i, 6);
		temp_scene.scene_orientation.kappa = Orientation_EOP(i, 7);


		// DEBUG
		cout << "EOP of ORIENTATION " << i << "\n\n";
		cout << "\tX:\t" << temp_scene.scene_orientation.X << endl;
		cout << "\tY:\t" << temp_scene.scene_orientation.Y << endl;
		cout << "\tZ:\t" << temp_scene.scene_orientation.Z << endl;
		cout << "\tOmega:\t" << temp_scene.scene_orientation.omega << endl;
		cout << "\tPhi:\t" << temp_scene.scene_orientation.phi << endl;
		cout << "\tKappa:\t" << temp_scene.scene_orientation.kappa << endl;


		// VISUALIZE
		//visualize_planes(planes);

		//Add to scene
		temp_scene.planes = planes;
		scenes.push_back(temp_scene);

		//DEBUG
		cout << "\n\nPlane equations for " << i << endl;
		for (int k = 0; k < planes.size(); k++)
		{
			cout << "\t" << planes[k].a1 << "\t" << planes[k].a2 << "\t" << planes[k].a3 << "\t" << planes[k].b << endl;
		}
		cout << "-------------------------- END -------------------------\n\n";
	}

	return scenes;

}

void plane_to_global(Plane &p1, Orientation O1)
{
	/*
		OMEGA: -180 to 180
		PHI: -180 to 180
		KAPPA: -180 to 180

	*/
	double del_omega, del_phi, del_kappa, plane_dist;
	Matrix3b3 R_del;
	RowVector3d global_translation, target_rot_vec, target_plane_vec, shiftdown;
	Vector4d full_transformed;

	target_plane_vec << p1.a1, p1.a2, p1.a3;

	Rotation_g2i(-1*O1.omega, -1*O1.phi, -1*O1.kappa, R_del);
	
	global_translation << O1.X, O1.Y, O1.Z;

	//Find new plane parameters
	target_rot_vec = target_plane_vec * R_del;

	//Find the distance. This is shiftdown coordinates

	plane_dist = -1*target_rot_vec * global_translation.transpose() + (p1.b);

	p1.a1 = target_rot_vec(0);
	p1.a2 = target_rot_vec(1);
	p1.a3 = target_rot_vec(2);
	p1.b = plane_dist;


}

Vector4d rotate_translate_plane(Matrix3b3 R, RowVector3d translation, Plane p1)
{
	Vector4d all_transformed;
	Matrix4b4 H;
	RowVector4d plane_equation;
	plane_equation << p1.a1, p1.a2, p1.a3, p1.b;
	H(0, 0) = R(0, 0); H(0, 1) = R(0, 1); H(0, 2) = R(0, 2); H(0, 3) = translation(0);
	H(1, 0) = R(1, 0); H(1, 1) = R(1, 1); H(1, 2) = R(1, 2); H(1, 3) = translation(1);
	H(2, 0) = R(2, 0); H(2, 1) = R(2, 1); H(2, 2) = R(2, 2); H(2, 3) = translation(2);
	H(3, 0) = 0;	   H(3, 1) = 0;		  H(3, 2) = 0;		 H(3, 3) = 1;

	all_transformed = H.inverse().transpose() * plane_equation.transpose();

	return all_transformed;

}

MatrixXd merge_data(MatrixXd IE_data, MatrixXd lidar_data, double time)
{
	MatrixXd output;
	output.resize(lidar_data.rows(), 11);
	double timestamp = 0;

	for (int i = 0; i < lidar_data.size(); i++)
	{
		output(i, 0) = i; //Point ID
		output(i, 2) = lidar_data(i, 12); // Horizontal angle
		output(i, 3) = lidar_data(i, 8); // Vertical angle
		output(i, 4) = lidar_data(i, 9); // Range

		for (int j = 0; j < IE_data.size(); j++)
		{
			if (IE_data(j, 0) == time) 
			{
				timestamp = IE_data(j, 0);
				output(i, 1) = timestamp;
				output(i, 5) = IE_data(j, 1); //X_GNSS
				output(i, 6) = IE_data(j, 2); //Y_GNSS
				output(i, 7) = IE_data(j, 3); //Z_GNSS
				output(i, 8) = IE_data(j, 7); //roll
				output(i, 9) = IE_data(j, 8); //pitch
				output(i, 10) = IE_data(j, 9); //yaw
			}
		}

	}

	return output;

}
void create_bundle_observations(vector<Scene> scenes, UniquePlanes unique, vector<RowVectorXd> &point_details, vector<RowVectorXd> &scene_details, vector<RowVectorXd> &plane_details)
{
	//This function makes the three Bundle Adjustment observation vectors

	RowVectorXd points_row(5), planes_row(4), scene_row(6);
	Scene temp_scene;
	double x, y, z;

	//Point Details
	//	This vector contains the details of each point in the adjustment. Each of these points are from a scene, 
	//	and lie on one of the unique planes
	//	| X | Y | Z | Unique Plane | Scene | 

	for (int i = 0; i < unique.mapping_vec.size(); i++)
	{
		temp_scene = scenes[unique.mapping_vec[i](1)];
		for (int j = 0; j < temp_scene.planes[unique.mapping_vec[i](2)].points_on_plane->points.size(); j++)
		{
			x = temp_scene.planes[unique.mapping_vec[i](2)].points_on_plane->points[j].x;
			y = temp_scene.planes[unique.mapping_vec[i](2)].points_on_plane->points[j].y;
			z = temp_scene.planes[unique.mapping_vec[i](2)].points_on_plane->points[j].z;
			points_row << x, y, z, unique.mapping_vec[i](0), unique.mapping_vec[i](1);
			point_details.push_back(points_row);
		}
	}

	//Plane Details
	//	This vector contains the details of each unique plane in the scenes.
	//	| A1 | A2 | A3 | B | 
	for (int i = 0; i < unique.unique_planes.size(); i++)
	{
		planes_row << unique.unique_planes[i].a1, unique.unique_planes[i].a2, unique.unique_planes[i].a3, unique.unique_planes[i].b;
		plane_details.push_back(planes_row);
	}

	//Scene Details
	//	This vector contains the details of each scene orientation details.
	//	| X | Y | Z | Omega | Phi | Kappa |
	for (int i = 0; i < scenes.size(); i++)
	{
		scene_row << scenes[i].scene_orientation.X, scenes[i].scene_orientation.Y, scenes[i].scene_orientation.Z, scenes[i].scene_orientation.omega, scenes[i].scene_orientation.phi, scenes[i].scene_orientation.kappa;
		scene_details.push_back(scene_row);
	}



}

vector<Scene> LoadDebugData()
{
	Scene temp_scene1, temp_scene2, temp_scene3, temp_scene4;
	vector<Scene> scenes;
	vector<char *> pcd_files1, pcd_files2, pcd_files3, pcd_files4;
	/*char * plane0_0 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O1_Planes\\Cloud_Plane_0.pcd";
	char * plane0_1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O1_Planes\\Cloud_Plane_1.pcd";
	char * plane0_2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O1_Planes\\Cloud_Plane_2.pcd";

	char * plane1_0 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O2_Planes\\Cloud_Plane_0.pcd";
	char * plane1_1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O2_Planes\\Cloud_Plane_1.pcd";
	char * plane1_2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O2_Planes\\Cloud_Plane_2.pcd";
	char * plane1_3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O2_Planes\\Cloud_Plane_2.pcd";

	char * plane2_0 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O3_Planes\\Cloud_Plane_0.pcd";
	char * plane2_1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O3_Planes\\Cloud_Plane_1.pcd";
	char * plane2_2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O3_Planes\\Cloud_Plane_2.pcd";

	char * plane3_0 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O4_Planes\\Cloud_Plane_0.pcd";
	char * plane3_1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O4_Planes\\Cloud_Plane_1.pcd";
	char * plane3_2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O4_Planes\\Cloud_Plane_2.pcd";

	char * plane_equations1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O1_Planes\\PlaneEquations.txt";
	char * plane_equations2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O2_Planes\\PlaneEquations.txt";
	char * plane_equations3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O3_Planes\\PlaneEquations.txt";
	char * plane_equations4 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O4_Planes\\PlaneEquations.txt";

	char * orientation1 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O1_Planes\\Orientation.txt";
	char * orientation2 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O2_Planes\\Orientation.txt";
	char * orientation3 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O3_Planes\\Orientation.txt";
	char * orientation4 = "C:\\Users\\Edmond\\Documents\\School\\Courses\\FifthYear\\ENGO500\\Data\\Crossiron\\Debug_INPUT\\O4_Planes\\Orientation.txt";*/
	

	char * plane0_0 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O1_Planes\\Cloud_Plane_0.pcd";
	char * plane0_1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O1_Planes\\Cloud_Plane_1.pcd";
	char * plane0_2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O1_Planes\\Cloud_Plane_2.pcd";
	
	char * plane1_0 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O2_Planes\\Cloud_Plane_0.pcd";
	char * plane1_1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O2_Planes\\Cloud_Plane_1.pcd";
	char * plane1_2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O2_Planes\\Cloud_Plane_2.pcd";
	char * plane1_3 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O2_Planes\\Cloud_Plane_2.pcd";
	
	char * plane2_0 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O3_Planes\\Cloud_Plane_0.pcd";
	char * plane2_1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O3_Planes\\Cloud_Plane_1.pcd";
	char * plane2_2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O3_Planes\\Cloud_Plane_2.pcd";
	
	char * plane3_0 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O4_Planes\\Cloud_Plane_0.pcd";
	char * plane3_1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O4_Planes\\Cloud_Plane_1.pcd";
	char * plane3_2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O4_Planes\\Cloud_Plane_2.pcd";
	
	char * plane_equations1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O1_Planes\\PlaneEquations.txt";
	char * plane_equations2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O2_Planes\\PlaneEquations.txt";
	char * plane_equations3 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O3_Planes\\PlaneEquations.txt";
	char * plane_equations4 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O4_Planes\\PlaneEquations.txt";
	
	char * orientation1 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O1_Planes\\Orientation.txt";
	char * orientation2 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O2_Planes\\Orientation.txt";
	char * orientation3 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O3_Planes\\Orientation.txt";
	char * orientation4 = "C:\\Users\\kiera.fulton2\\Desktop\\ENGO500\\CmakeTests\\Test2\\Build\\Debug_INPUT\\O4_Planes\\Orientation.txt";


	pcd_files1.push_back(plane0_0);
	pcd_files1.push_back(plane0_1);
	pcd_files1.push_back(plane0_2);


	pcd_files2.push_back(plane1_0);
	pcd_files2.push_back(plane1_1);
	pcd_files2.push_back(plane1_2);


	pcd_files3.push_back(plane2_0);
	pcd_files3.push_back(plane2_1);
	pcd_files3.push_back(plane2_2);

	pcd_files4.push_back(plane3_0);
	pcd_files4.push_back(plane3_1);
	pcd_files4.push_back(plane3_2);


	temp_scene1.planes = get_debug_planes(plane_equations1);
	temp_scene1.scene_orientation = get_debug_orientation(orientation1);

	PointCloudXYZptr plane_cloud1_1(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud1_2(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud1_3(new PointCloudXYZ);
	Read_Lidar_points(pcd_files1[0], plane_cloud1_1);
	Read_Lidar_points(pcd_files1[1], plane_cloud1_2);
	Read_Lidar_points(pcd_files1[2], plane_cloud1_3);
	temp_scene1.planes[0].points_on_plane = plane_cloud1_1;
	temp_scene1.planes[1].points_on_plane = plane_cloud1_2;
	temp_scene1.planes[2].points_on_plane = plane_cloud1_3;
	scenes.push_back(temp_scene1);
	//------------------------------------------------------------------
	temp_scene2.planes = get_debug_planes(plane_equations2);
	temp_scene2.scene_orientation = get_debug_orientation(orientation2);

	PointCloudXYZptr plane_cloud2_1(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud2_2(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud2_3(new PointCloudXYZ);
	Read_Lidar_points(pcd_files2[0], plane_cloud2_1);
	Read_Lidar_points(pcd_files2[1], plane_cloud2_2);
	Read_Lidar_points(pcd_files2[2], plane_cloud2_3);
	temp_scene2.planes[0].points_on_plane = plane_cloud2_1;
	temp_scene2.planes[1].points_on_plane = plane_cloud2_2;
	temp_scene2.planes[2].points_on_plane = plane_cloud2_3;
	scenes.push_back(temp_scene2);
	//----------------------------------------------------------------------
	temp_scene3.planes = get_debug_planes(plane_equations3);
	temp_scene3.scene_orientation = get_debug_orientation(orientation3);

	PointCloudXYZptr plane_cloud3_1(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud3_2(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud3_3(new PointCloudXYZ);
	Read_Lidar_points(pcd_files3[0], plane_cloud3_1);
	Read_Lidar_points(pcd_files3[1], plane_cloud3_2);
	Read_Lidar_points(pcd_files3[2], plane_cloud3_3);
	temp_scene3.planes[0].points_on_plane = plane_cloud3_1;
	temp_scene3.planes[1].points_on_plane = plane_cloud3_2;
	temp_scene3.planes[2].points_on_plane = plane_cloud3_3;
	scenes.push_back(temp_scene3);

	//-------------------------------------------------------------------
	temp_scene4.planes = get_debug_planes(plane_equations4);
	temp_scene4.scene_orientation = get_debug_orientation(orientation4);

	PointCloudXYZptr plane_cloud4_1(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud4_2(new PointCloudXYZ);
	PointCloudXYZptr plane_cloud4_3(new PointCloudXYZ);
	Read_Lidar_points(pcd_files4[0], plane_cloud4_1);
	Read_Lidar_points(pcd_files4[1], plane_cloud4_2);
	Read_Lidar_points(pcd_files4[2], plane_cloud4_3);
	temp_scene4.planes[0].points_on_plane = plane_cloud4_1;
	temp_scene4.planes[1].points_on_plane = plane_cloud4_2;
	temp_scene4.planes[2].points_on_plane = plane_cloud4_3;
	scenes.push_back(temp_scene4);

	return scenes;
}

vector<Plane> get_debug_planes(char *filename)
{
	Plane temp_plane;
	vector<Plane> planes;
	MatrixXd inmat;
	Read_Mat(filename, inmat);

	for (int i = 0; i < inmat.rows(); i++)
	{
		temp_plane.a1 = inmat(i, 0);
		temp_plane.a2 = inmat(i, 1);
		temp_plane.a3 = inmat(i, 2);
		temp_plane.b = inmat(i, 3);
		planes.push_back(temp_plane);
	}

	return planes;

}


Orientation get_debug_orientation(char *filename)
{
	Orientation temp_o;
	MatrixXd inmat;
	Read_Mat(filename, inmat);

	temp_o.X = inmat(0,0);
	temp_o.Y = inmat(0, 1);
	temp_o.Z = inmat(0, 2);
	temp_o.omega = inmat(0, 3);
	temp_o.phi = inmat(0, 4);
	temp_o.kappa = inmat(0, 5);

	return temp_o;

}


double check_plane_az(Orientation base_O, Orientation target_O, Plane plane_base, Plane plane_target)
{
	double base_az = get_plane_az(plane_base); //local azimuth (rad)
	double target_az = get_plane_az(plane_target); //local azimuth (rad)
	
	//change to degree
	base_az = base_az*RAD2DEG;
	target_az = target_az*RAD2DEG;


	base_az = base_az + base_O.kappa; //Change to global
	target_az = target_az + target_O.kappa; //Change to global

	if (base_az > 360) { base_az = base_az - 360; }
	if (target_az > 360) { target_az = target_az - 360; }

	return abs(target_az - base_az);


}


double get_plane_az(Plane test_plane)
{
	// | The azimuth of a point is defined as the clockwise angle from north. In this case the y axis

	double temp_az, sumx, sumy, meanx, meany;
	double az_accumulate = 0;
	sumx = 0;
	sumy = 0;
	for (int i = 0; i < test_plane.points_on_plane->size(); i++)
	{

		sumx = sumx + test_plane.points_on_plane->points[i].x;
		sumy = sumy + test_plane.points_on_plane->points[i].y;

		//temp_az = atan2(test_plane.points_on_plane->points[i].x, test_plane.points_on_plane->points[i].y);
		//if (temp_az < 0) { temp_az = (360*1/RAD2DEG) + temp_az; }
		//az_accumulate = az_accumulate + temp_az;
		////clog << "The azimuth of point " << i << ": x " << test_plane.points_on_plane->points[i].x << " y: " << test_plane.points_on_plane->points[i].y << "  is: " << temp_az*RAD2DEG << endl;
	}
	meanx = sumx / test_plane.points_on_plane->size();
	meany = sumy / test_plane.points_on_plane->size();

	//combined_data
	temp_az = atan2(meanx , meany);

	if (temp_az < 0) { temp_az = (360 * 1 / RAD2DEG) + temp_az; }

	return temp_az;
}

void get_hour_day(double GPS_time, double *Hour, int *Day)
{

	if (GPS_time >= 0 && GPS_time < 86400)
	{
		*Day = 0;
	}
	else if (GPS_time >= 86400 && GPS_time < 172800)
	{
		*Day = 1;
	}
	else if (GPS_time >= 172800 && GPS_time < 259200)
	{
		*Day = 2;
	}
	else if (GPS_time >= 259200 && GPS_time < 345600)
	{
		*Day = 3;
	}
	else if (GPS_time >= 345600 && GPS_time < 432000)
	{
		*Day = 4;
	}
	else if (GPS_time >= 432000 && GPS_time < 518400)
	{
		*Day = 5;
	}
	else if (GPS_time >= 518400)
	{
		*Day = 6;
	}
	else
	{
		cout << "Problem with first GPS time. " << endl;
	}

	double temp = (GPS_time - (86400 * *Day)) / 3600;
	*Hour = floor(temp);
}

double round_time(double time)
{
	double whole, frac;
	whole = floor(time);
	frac = time - whole;

	//Round to two decimal places
	double two_dec = floor(time * 100) / 100;

	//Round to nearest quarter integer
	double rounded = floor((time * 4) + 0.25) / 4;
	return rounded;
}