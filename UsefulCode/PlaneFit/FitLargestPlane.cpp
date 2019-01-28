
// SOURCE: https://stackoverflow.com/questions/46826720/pclransac-segmentation-get-all-planes-in-cloud
// USEFUL: 
// 			https://stackoverflow.com/questions/23814632/ransac-plane-fitting-coefficients
//			http://pointclouds.org/documentation/tutorials/planar_segmentation.php			

int
main(int argc, char** argv)
{
    pcl::visualization::CloudViewer viewer("viewer1");

    pcl::PCLPointCloud2::Ptr cloud_blob(new pcl::PCLPointCloud2), cloud_filtered_blob(new pcl::PCLPointCloud2);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>), cloud_p(new pcl::PointCloud<pcl::PointXYZ>), cloud_f(new pcl::PointCloud<pcl::PointXYZ>);

    // Fill in the cloud data
    pcl::PCDReader reader;
    reader.read("clouds/table.pcd", *cloud_blob);

    // Create the filtering object: downsample the dataset using a leaf size of 1cm
    pcl::VoxelGrid<pcl::PCLPointCloud2> sor;
    sor.setInputCloud(cloud_blob);
    sor.setLeafSize(0.01f, 0.01f, 0.01f);
    sor.filter(*cloud_filtered_blob);

    // Convert to the templated PointCloud
    pcl::fromPCLPointCloud2(*cloud_filtered_blob, *cloud_filtered);

    std::cerr << "PointCloud after filtering: " << cloud_filtered->width * cloud_filtered->height << " data points." << std::endl;

    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
    // Create the segmentation object
    pcl::SACSegmentation<pcl::PointXYZ> seg;
    // Optional
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_PLANE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setMaxIterations(1000);
    seg.setDistanceThreshold(0.01);

    // Create the filtering object
    pcl::ExtractIndices<pcl::PointXYZ> extract;

    int i = 0, nr_points = (int)cloud_filtered->points.size();
    // While 30% of the original cloud is still there
    while (cloud_filtered->points.size() > 0.3 * nr_points)
    {
        // Segment the largest planar component from the remaining cloud
        seg.setInputCloud(cloud_filtered);
        pcl::ScopeTime scopeTime("Test loop");
        {
            seg.segment(*inliers, *coefficients);
        }
        if (inliers->indices.size() == 0)
        {
            std::cerr << "Could not estimate a planar model for the given dataset." << std::endl;
            break;
        }

        // Extract the inliers
        extract.setInputCloud(cloud_filtered);
        extract.setIndices(inliers);
        extract.setNegative(false);
        extract.filter(*cloud_p);
        std::cerr << "PointCloud representing the planar component: " << cloud_p->width * cloud_p->height << " data points." << std::endl;

        }

    viewer.showCloud(cloud_p, "viewer1");
    while (!viewer.wasStopped()) {}

    return (0);
}




	clog << "\n-------------------------STEP 2: upsample data -------------------------------------------------------\n";

	/*

	// Create K-D Tree
	//pcl::search::KdTree<pcl::PointXYZI>::Ptr kd_tree(new pcl::search::KdTree<pcl::PointXYZI>);

	// Create point normals
	//pcl::PointCloud<pcl::PointNormal> normals;

	// Create a KD-Tree
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);

	// Output has the PointNormal type in order to store the normals calculated by MLS
	pcl::PointCloud<pcl::PointNormal> mls_points;

	// Init object (second point type is for the normals, even if unused)
	pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
	
	mls.setComputeNormals(true);
	
	// Set parameters
	mls.setInputCloud(Novatel_cloud);
	mls.setPolynomialOrder(2);
	mls.setSearchMethod(tree);
	mls.setSearchRadius(0.03);
	mls.setUpsamplingMethod(pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal>::UpsamplingMethod::SAMPLE_LOCAL_PLANE);
	mls.setUpsamplingRadius(0.005);
	mls.setUpsamplingStepSize(0.005);
	
	// Reconstruct
	mls.process(mls_points);
	
	//Place back in point cloud
	pcl::copyPointCloud(mls_points, *Novatel_cloud);

	pcl::PCDWriter writer;
	writer.write("Upsampled_PointCloud.pcd", *Novatel_cloud);

	*/