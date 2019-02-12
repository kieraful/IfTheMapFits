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

		std::clog << "File read successful\n";
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
		0, cos(Omega), sin(Omega),
		0, -sin(Omega), cos(Omega);

	Mf0 << cos(Phi), 0, -sin(Phi),
		0, 1, 0,
		sin(Phi), 0, cos(Phi);

	Mk0 << cos(Kappa), sin(Kappa), 0,
		-sin(Kappa), cos(Kappa), 0,
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

	double percent_cloud = 50;

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
	Eigen::Vector3f search_axis; // Axis to search for planes PERPENDICULAR to
	double plane_buffer = 22.5*PI / 180;//Degree offset from search plane to allow
	double size, remaining_p, remaining_pts;

	// Fill in the cloud data
	pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
	

	// Create the segmentation object
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	// Optional
	//Set search axis
	//search_axis << 0, 1, 0; //y axis
	//seg.setOptimizeCoefficients(true);
	//seg.setModelType(pcl::SACMODEL_PLANE);
	//seg.setMethodType(pcl::SAC_RANSAC);
	//seg.setMaxIterations(500);
	//seg.setAxis(search_axis);
	//seg.setEpsAngle(plane_buffer);
	//seg.setDistanceThreshold(0.01);

	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PARALLEL_PLANE); //only want points perpendicular to a given axis
	seg.setMaxIterations(1000);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setDistanceThreshold(0.15); // keep points within 0.10 m of the plane
	Eigen::Vector3f axis = Eigen::Vector3f(0.0, 0.0, 1.0); //x axis
	seg.setAxis(axis);
	seg.setEpsAngle(20.0f * (PI / 180.0f)); // plane can be within 30 degrees of X-Z plane
	

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

		//cout << "\n all pts: " << nr_points << " remaining percent: " << remaining_p << " Planed pts: " << remaining_pts;
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
		//cerr << "The plane " << n_planes <<" has " << cloud_p->width * cloud_p->height << " data points." << endl;
		//cerr << "\tThe plane has coefficiants: a= " << temp_plane.a1 << " b= " << temp_plane.a2 << " c= " << temp_plane.a3 << " \n";

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

void save_planes(vector<Plane> planes)
{

	pcl::PCDWriter writer; //writer object for point clouds

	for (int i = 0; i < planes.size(); i++)
	{
		std::stringstream ss;
		ss << "Cloud_Plane_" << i << ".pcd";
		writer.write<pcl::PointXYZ>(ss.str(), *planes[i].points_on_plane, false);

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
	output.resize(num_points, 4);


	double timestamp = 0.0;
	double vert_angle = 0.0;
	double horiz_angle = 0.0;
	double range = 0.0;

	//Iterating for each epoch
	for (int i = 0; i < num_points; i++)
	{
		timestamp = data(i, 0);

		//Populate matrices with applicable values
		r_b_geo << data(i, 1),
			data(i, 2),
			data(i, 3);

		Rotation_g2i(data(i, 7), data(i, 8), data(i, 9), R_b_geo); // Inputs: roll, pitch, azimuth, output matrix

		r_lidar_b << boresight_LA(0, 0),
			boresight_LA(0, 1),
			boresight_LA(0, 2); // Assumes boresight_LA read in is a 1x3 matrix 

		R_lidar_b << boresight_angles(0, 0),
			boresight_angles(0, 1),
			boresight_angles(0, 2); // Assumes boresight_angles is a 1x3 matrix -> read these into Rotation_g2i matrix

		Rotation_g2i(boresight_angles(0, 0), boresight_angles(0, 1), boresight_angles(0, 2), R_lidar_b); // Assumes boresight_angles is a 1x3 matrix
			

		vert_angle = data(i, 13);
		horiz_angle = data(i, 14);
		range = data(i, 15);

		r_p_lidar << range * cos(vert_angle) * sin(horiz_angle),
			range * cos(vert_angle) * cos(horiz_angle),
			range*sin(vert_angle);

		r_p_geo = r_b_geo + (R_b_geo * r_lidar_b) + (R_b_geo * R_lidar_b * r_p_lidar);

		// push r_p_geo on to output with timestamp
		output(i, 0) = timestamp;
		output(i, 1) = r_p_geo(0, 0);
		output(i, 2) = r_p_geo(1, 0);
		output(i, 3) = r_p_geo(2, 0);

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
	RowVector3d target_plane_vec, base_plane_vec, target_rot_vec, mapping_temp;
	vector<int> candidates;
	double thresh_orientation = 0.1;
	double best_dist, dot_prod, dist_temp;
	int best_plane;
	Plane temp_plane;
	Orientation temp_orientation;

	//Add base scene to unique planes, as there is nothing to match yet
	for (int i = 0; i < scenes[0].planes.size(); i++) {

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
			// Target plane to matrix
			target_plane_vec << scenes[i].planes[j].a1, scenes[i].planes[j].a2, scenes[i].planes[j].a3;

			candidates.clear();
			for (int k = 0; k < unique.unique_planes.size();k++) // for each unique plane k
			{
				// Make rotation matrix from scene i to unique plane
				del_omega = scenes[i].scene_orientation.omega - unique.reference_orientations[k].omega;
				del_phi = scenes[i].scene_orientation.phi - unique.reference_orientations[k].phi;
				del_kappa = scenes[i].scene_orientation.kappa - unique.reference_orientations[k].kappa;
				Rotation_g2i(del_omega, del_kappa, del_phi, R_del);

				// Rotate plane to match base
				target_rot_vec = target_plane_vec * R_del;

				// Base plane to matrix
				base_plane_vec << unique.unique_planes[k].a1, unique.unique_planes[k].a2, unique.unique_planes[k].a3;
				//dot product

				 //------ Debug --------------

				//cout << "\n\n \t\tDEBUGGING\t\n\n";
				//cout << "The DOT product: " << base_plane_vec.dot(target_plane_vec) << endl;
				//---------------------------


				dot_prod = base_plane_vec.dot(target_plane_vec);
				if (1 - dot_prod < thresh_orientation)
				{
					candidates.push_back(k); // place candidate unique plane index in vector
				}


			}

			// If there is more than 1 candidate, check the relative distances

			// reset distance threshold
			best_dist = 10; // Planes should definitely not be more than 2 meters away from each other
			best_plane = -1;
			// For each candidate plane, find closest matching plane in base (Euclidian distance)
			if (candidates.size() > 1)
			{
				for (int m = 0; m < candidates.size(); m++)
				{
					//dist_temp = abs(scenes[i].planes[k].b - unique.unique_planes[candidates[m]].b);

					//clog << "Checking distance!";



					dist_temp = check_plane_dists(unique.reference_orientations[candidates[m]], scenes[i].scene_orientation, unique.unique_planes[candidates[m]], scenes[i].planes[j]);
					//clog << dist_temp << endl;

					//Debug
					//cout << "\nDistance : " << dist_temp;
					//cout << "\tm: " << m;
					//cout << "\tCandidate[m]" << candidates[m];
					//--------------

					if (dist_temp < best_dist)
					{
						//This is the best plane so far
						best_dist = dist_temp;
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
				//clog << "\nFound matching plane for target scene " << i << " plane " << k << ", in base " << j << " plane " << best_plane << endl;
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

void print_vector(vector<int> print_vector)
{
	for (int i = 0; i < print_vector.size(); i++)
	{

		cout << "\t" << print_vector[i] << endl;

	}

}


double check_plane_dists(Orientation orient_base, Orientation orient_target, Plane plane_base, Plane plane_target)
{
	// Checks the distance between two planes to see if they are in fact the same plane
	// Done by finding circle intersection, and if it satisfies the base plane equation. 

	bool result = false;
	double d_base, d_target; 
	double x1, x2, y1, y2, r1, r2, d, a, h, px, py, int1x, int1y, int2x, int2y, check1, check2;


	d_base = abs(plane_base.b); // Distance from origin in base scene
	d_target = abs(plane_target.b); // Distance from origin in target scene

	//find intersection points
	x1 = orient_base.X;
	y1 = orient_base.Y;
	x2 = orient_target.X;
	y2 = orient_target.X;

	r1 = abs(plane_base.b);
	r2 = abs(plane_target.b);

	d = sqrt(pow((x1 - x2),2) + pow(y1 - y2,2));

	if (d > (r1 + r2))
	{
		//Does not intersect. Not the right plane. 
		check1 = 100;
		check2 = 100;

	}
	else if (d + min(r1, r2) <= max(r1, r2))
	{

		// One is in other
		//TODO: Find closest point

		check1 = abs(r1 - r2 - d);
		check2 = 100;



	}

	else
	{
		a = (pow(r1,2) - pow(r2, 2) + pow(d,2)) / (2 * d);
		h = sqrt(pow(r1, 2) - pow(a, 2));

		px = x1 + (a*(x2 - x1)) / d;
		py = y1 + (a*(y2 - y1)) / d;

		int1x = px + (h*(y2 - y1)) / d;
		int1y = py - (h*(x2 - x1)) / d;

		int2x = px - (h*(y2 - y1)) / d;
		int2y = py + (h*(x2 - x1)) / d;
		

		// Now we have the intersection. Check if either one lies on the plane. 
		check1 = abs(plane_base.a1*int1x + plane_base.a2*int1y + plane_base.a3*orient_target.Z - plane_base.b);
		check2 = abs(plane_base.a1*int2x + plane_base.a2*int2x + plane_base.a3*orient_target.Z - plane_base.b);

	}


	return min(check1, check2);
}