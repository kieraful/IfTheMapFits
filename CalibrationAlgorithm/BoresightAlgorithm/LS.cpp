#include "LS.h"

//Receives the name of a file "FileName" containing a numerical matrix, and read the matrix data into variable "m"
//void Read_Mat(char *FileName, MatrixXd& m) {
//
//	m.resize(0, 0);
//
//	ifstream matfile;
//	matfile.open(FileName, ios::in); //to open the file and start reading from the beginning
//	if (matfile.fail()) //check if the file is opened successfully
//	{
//		cout << "\nThere was a problem reading the following file: " << endl << FileName << endl;
//		//exit (EXIT_FAILURE);
//		return;
//	}
//	else {
//		cout << "\nFile read correctly. Continuing.\n";
//	}
//
//	char* readlinechr = new char[MaxMatSize];
//	vector<double> v_all;
//	int nrow = 0;
//
//	while (matfile.getline(readlinechr, MaxMatSize, '\n')) {
//		nrow++;
//		int stln = strlen(readlinechr);
//		char* readlinestr = new char[stln + 1];
//		for (int i = 0; i<stln; i++)
//		{
//			readlinestr[i] = readlinechr[i];
//		}
//
//		readlinestr[stln] = '\0';
//
//		stringstream rowstream(readlinestr);
//		double value;
//		while (!rowstream.eof()) {
//			rowstream >> value;
//			v_all.push_back(value);
//		}
//	}
//	matfile.close();
//
//	int ncol = v_all.size() / nrow;
//	m.resize(nrow, ncol);
//
//	for (int i = 0; i<nrow; i++) {
//		for (int j = 0; j<ncol; j++) {
//			m(i, j) = v_all.at(i*ncol + j);
//		}
//	}
//
//	return;
//}

//Function to fill a matrix with info from vector of vectors
void vec2mat(vector<RowVectorXd>& vec, MatrixXd& mat, int cols)
{
	mat = MatrixXd(vec.size(), cols);
	for (int i = 0; i<vec.size(); i++)
	{
		for (int j = 0; j<cols; j++)
		{
			mat(i, j) = vec[i][j];
		}
	}

	return;
}



//Function to compute A
void computeAandw(MatrixXd& A_full, MatrixXd& A, MatrixXd& H, MatrixXd& w_full, MatrixXd& w, MatrixXd& V, int u, int numPlanes, int numScans, int numLidPts, MatrixXd& bs_params, MatrixXd& plane_details, MatrixXd& scene_details, MatrixXd& point_details, MatrixXd& GNSS_INS_data)
{
	//bs_params: initial approx boresight parameters for entire adjustment. 6 by 1.
	//x_Sjb y_Sjb z_Sjb w_Sb phi_Sb K_Sb

	//plane_details: normals and distance for each plane in vector of vectors. numPlanes by 4.
	//n_xpg, n_ypg, n_zpg, d_p
	//a1, a2, a3, b

	//scene_details: 6 position and orientation parameters from lidar to ground for each scan in vector of vectors. numScans by 6.
	//x_bjg y_bjg z_bjg w_bjg phi_bjg K_bjg
	//X, Y, Z, omega, phi, kappa;

	//point_details format: all point coordinates OBSERVED with lidar in vector of vectors. numLidPts by 3.
	//The same actual point can be observed in diff scenes and have diff x,y,z values.
	//x_sj y_sj z_sj
	//x y z


	MatrixXd A1(0, 0);
	MatrixXd A2(0, 0);
	MatrixXd A3(0, 0);

	MatrixXd w1(0, 0);
	MatrixXd w2(0, 0);
	MatrixXd w3(0, 0);


	//Loop for scans
	for (int i = 0; i<numScans; i++) {
		A1.conservativeResize(A1.rows() + 6, u);
		A1.block(A1.rows() - 6, 0, 6, u) = computeAScan(u, numPlanes, i);

		w1.conservativeResize(w1.rows() + 6, 1);
		w1.block(w1.rows() - 6, 0, 6, 1) = computewScan(scene_details(i, 0), scene_details(i, 1), scene_details(i, 2),
			scene_details(i, 3), scene_details(i, 4), scene_details(i, 5),
			GNSS_INS_data(i, 0), GNSS_INS_data(i, 1), GNSS_INS_data(i, 2),
			GNSS_INS_data(i, 3), GNSS_INS_data(i, 4), GNSS_INS_data(i, 5));
	}

	//cout << endl << A1 << endl << endl;


	for (int j = 0; j<numPlanes; j++) {
		A2.conservativeResize(A2.rows() + 1, u);
		A2.block(A2.rows() - 1, 0, 1, u) = computeAPlane(u, j, plane_details(j, 0), plane_details(j, 1), plane_details(j, 2));

		w2.conservativeResize(w2.rows() + 1, 1);
		w2.block(w2.rows() - 1, 0, 1, 1) = computewPlane(plane_details(j, 0), plane_details(j, 1), plane_details(j, 2));
	}
	//cout << A2 << endl << endl;


	int planeID = 0;
	int scanID = 0;
	for (int k = 0; k<numLidPts; k++) {

		A3.conservativeResize(A3.rows() + 1, u);
		w3.conservativeResize(w3.rows() + 1, 1);

		planeID = point_details(k, 3);
		scanID = point_details(k, 4);

		A3.block(A3.rows() - 1, 0, 1, u) = computeAPt(u, numPlanes, planeID, scanID,
			bs_params(0, 0), bs_params(1, 0), bs_params(2, 0),
			bs_params(3, 0), bs_params(4, 0), bs_params(5, 0),
			plane_details(planeID, 0), plane_details(planeID, 1), plane_details(planeID, 2),
			scene_details(scanID, 0), scene_details(scanID, 1), scene_details(scanID, 2),
			scene_details(scanID, 3), scene_details(scanID, 4), scene_details(scanID, 5),
			point_details(k, 0), point_details(k, 1), point_details(k, 2));

		w3.block(w3.rows() - 1, 0, 1, 1) = computewPt(bs_params(0, 0), bs_params(1, 0), bs_params(2, 0),
			bs_params(3, 0), bs_params(4, 0), bs_params(5, 0),
			plane_details(planeID, 0), plane_details(planeID, 1), plane_details(planeID, 2), plane_details(planeID, 3),
			scene_details(scanID, 0), scene_details(scanID, 1), scene_details(scanID, 2),
			scene_details(scanID, 3), scene_details(scanID, 4), scene_details(scanID, 5),
			point_details(k, 0), point_details(k, 1), point_details(k, 2));
	}

	//for debugging
	//cout << i << "," << j << "," << k << "\n";


	//Combine scan rows, plane rows, and column rows into the full A matrix
	//cout << A_full.rows() << "," << A_full.cols() << "\n";
	//cout << A1.rows() << "," << A1.cols() << "\n";
	//cout << A2.rows() << "," << A2.cols() << "\n";
	//cout << A3.rows() << "," << A3.cols() << "\n";
	A_full.block(0, 0, A1.rows(), u) = A1;
	A_full.block(numScans * 6, 0, A2.rows(), u) = A2;
	A_full.block(numScans * 6 + numPlanes, 0, A3.rows(), u) = A3;

	A.block(0, 0, A1.rows(), u) = A1;
	A.block(numScans * 6, 0, A3.rows(), u) = A3;

	H = A2;

	A1.resize(0, 0);
	A2.resize(0, 0);
	A3.resize(0, 0);


	w_full.block(0, 0, w1.rows(), 1) = w1;
	w_full.block(numScans * 6, 0, w2.rows(), 1) = w2;
	w_full.block(numScans * 6 + numPlanes, 0, w3.rows(), 1) = w3;

	w.block(0, 0, w1.rows(), 1) = w1;
	w.block(numScans * 6, 0, w3.rows(), 1) = w3;

	V = w2;


	return;
}



//Function to compute 6 rows of A for a scan
//scanNum is the number of the scan that matches the scan equations
MatrixXd computeAScan(int u, int numPlanes, int scanNum)
{
	MatrixXd A_scan = MatrixXd::Zero(6, u);

	int i = 6 + 4 * numPlanes + 6 * scanNum;
	int k = 0;

	for (int j = i; k<6; j++) {
		A_scan(k, j) = 1;
		k++;
	}

	return A_scan;
}



//Function to compute a row of A for a plane
//planeNum is the number of the plane that matches the plane equation
MatrixXd computeAPlane(int u, int planeNum, double n_xpg, double n_ypg, double n_zpg)
{
	MatrixXd A_plane = MatrixXd::Zero(1, u);

	int i = 6 + 4 * planeNum;

	A_plane(0, i) = 2 * n_xpg;
	A_plane(0, i + 1) = 2 * n_ypg;
	A_plane(0, i + 2) = 2 * n_zpg;

	return A_plane;
}


//Function to compute elements of a rotation matrix, image to object rotation
MatrixXd RotMatElements(double w, double phi, double K)
{
	MatrixXd Rotw(3, 3);
	MatrixXd Rotphi(3, 3);
	MatrixXd RotK(3, 3);

	Rotw << 1, 0, 0,
		0, cosd(w), sind(w),
		0, -sind(w), cosd(w);

	Rotphi << cosd(phi), 0, -sind(phi),
		0, 1, 0,
		sind(phi), 0, cosd(phi);

	RotK << cosd(K), sind(K), 0,
		-sind(K), cosd(K), 0,
		0, 0, 1;

	MatrixXd RotMat = RotK*Rotphi*Rotw;
	RotMat.transposeInPlace();

	return RotMat;
}



//Function to compute the derivatives of a point equation wrt rotation matrix elements for Rbjg
MatrixXd PtEqnWrtRotbjg(double x_Sjb, double y_Sjb, double z_Sjb, double w_Sb, double phi_Sb, double K_Sb, double x_sj, double y_sj, double z_sj, double n_xpg, double n_ypg, double n_zpg)
{
	MatrixXd Dbjg1(9, 1);

	MatrixXd RotMat = RotMatElements(w_Sb, phi_Sb, K_Sb);

	double RSb_11 = RotMat(0, 0);
	double RSb_12 = RotMat(0, 1);
	double RSb_13 = RotMat(0, 2);
	double RSb_21 = RotMat(1, 0);
	double RSb_22 = RotMat(1, 1);
	double RSb_23 = RotMat(1, 2);
	double RSb_31 = RotMat(2, 0);
	double RSb_32 = RotMat(2, 1);
	double RSb_33 = RotMat(2, 2);

	Dbjg1(0, 0) = n_xpg*(x_Sjb + RSb_11*x_sj + RSb_12*y_sj + RSb_13*z_sj);
	Dbjg1(1, 0) = n_xpg*(y_Sjb + RSb_21*x_sj + RSb_22*y_sj + RSb_23*z_sj);
	Dbjg1(2, 0) = n_xpg*(z_Sjb + RSb_31*x_sj + RSb_32*y_sj + RSb_33*z_sj);
	Dbjg1(3, 0) = n_ypg*(x_Sjb + RSb_11*x_sj + RSb_12*y_sj + RSb_13*z_sj);
	Dbjg1(4, 0) = n_ypg*(y_Sjb + RSb_21*x_sj + RSb_22*y_sj + RSb_23*z_sj);
	Dbjg1(5, 0) = n_ypg*(z_Sjb + RSb_31*x_sj + RSb_32*y_sj + RSb_33*z_sj);
	Dbjg1(6, 0) = n_zpg*(x_Sjb + RSb_11*x_sj + RSb_12*y_sj + RSb_13*z_sj);
	Dbjg1(7, 0) = n_zpg*(y_Sjb + RSb_21*x_sj + RSb_22*y_sj + RSb_23*z_sj);
	Dbjg1(8, 0) = n_zpg*(z_Sjb + RSb_31*x_sj + RSb_32*y_sj + RSb_33*z_sj);

	return Dbjg1;
}



//Function to compute the derivatives of a point equation wrt rotation matrix elements for RSb
MatrixXd PtEqnWrtRotSb(double w_bjg, double phi_bjg, double K_bjg, double x_sj, double y_sj, double z_sj, double n_xpg, double n_ypg, double n_zpg)
{
	MatrixXd DSb1(9, 1);

	MatrixXd RotMat = RotMatElements(w_bjg, phi_bjg, K_bjg);

	double Rbjg_11 = RotMat(0, 0);
	double Rbjg_12 = RotMat(0, 1);
	double Rbjg_13 = RotMat(0, 2);
	double Rbjg_21 = RotMat(1, 0);
	double Rbjg_22 = RotMat(1, 1);
	double Rbjg_23 = RotMat(1, 2);
	double Rbjg_31 = RotMat(2, 0);
	double Rbjg_32 = RotMat(2, 1);
	double Rbjg_33 = RotMat(2, 2);

	DSb1(0, 0) = x_sj*(Rbjg_11*n_xpg + Rbjg_21*n_ypg + Rbjg_31*n_zpg);
	DSb1(1, 0) = y_sj*(Rbjg_11*n_xpg + Rbjg_21*n_ypg + Rbjg_31*n_zpg);
	DSb1(2, 0) = z_sj*(Rbjg_11*n_xpg + Rbjg_21*n_ypg + Rbjg_31*n_zpg);
	DSb1(3, 0) = x_sj*(Rbjg_12*n_xpg + Rbjg_22*n_ypg + Rbjg_32*n_zpg);
	DSb1(4, 0) = y_sj*(Rbjg_12*n_xpg + Rbjg_22*n_ypg + Rbjg_32*n_zpg);
	DSb1(5, 0) = z_sj*(Rbjg_12*n_xpg + Rbjg_22*n_ypg + Rbjg_32*n_zpg);
	DSb1(6, 0) = x_sj*(Rbjg_13*n_xpg + Rbjg_23*n_ypg + Rbjg_33*n_zpg);
	DSb1(7, 0) = y_sj*(Rbjg_13*n_xpg + Rbjg_23*n_ypg + Rbjg_33*n_zpg);
	DSb1(8, 0) = z_sj*(Rbjg_13*n_xpg + Rbjg_23*n_ypg + Rbjg_33*n_zpg);

	return DSb1;
}



//Function to compute derivatives of rotation matrix elements wrt rotation angles
MatrixXd RotWrtAngles(double w, double phi, double K)
{
	MatrixXd D = MatrixXd::Zero(9, 3);

	D(0, 1) = -cosd(K)*sind(phi);
	D(0, 2) = -sind(K)*cosd(phi);
	//[0, -cosd(K)*sind(phi), -sind(K)*cosd(phi)]

	D(1, 0) = cosd(K)*cosd(w)*sind(phi) - sind(K)*sind(w);
	D(1, 1) = cosd(K)*cosd(phi)*sind(w);
	D(1, 2) = cosd(K)*cosd(w) - sind(K)*sind(phi)*sind(w);
	//[cosd(K)*cosd(w)*sind(phi) - sind(K)*sind(w), cosd(K)*cosd(phi)*sind(w), cosd(K)*cosd(w) - sind(K)*sind(phi)*sind(w)]

	D(2, 0) = sind(K)*cosd(w) + cosd(K)*sind(phi)*sind(w);
	D(2, 1) = -cosd(K)*cosd(phi)*cosd(w);
	D(2, 2) = cosd(K)*sind(w) + sind(K)*cosd(w)*sind(phi);
	//[sind(K)*cosd(w) + cosd(K)*sind(phi)*sind(w), -cosd(K)*cosd(phi)*cosd(w), cosd(K)*sind(w) + sind(K)*cosd(w)*sind(phi)]

	D(3, 1) = sind(K)*sind(phi);
	D(3, 2) = -cosd(K)*cosd(phi);
	//[0, sind(K)*sind(phi), -cosd(K)*cosd(phi)]

	D(4, 0) = -cosd(K)*sind(w) - sind(K)*cosd(w)*sind(phi);
	D(4, 1) = -sind(K)*cosd(phi)*sind(w);
	D(4, 2) = -sind(K)*cosd(w) - cosd(K)*sind(phi)*sind(w);
	//[ - cosd(K)*sind(w) - sind(K)*cosd(w)*sind(phi), -sind(K)*cosd(phi)*sind(w), - sind(K)*cosd(w) - cosd(K)*sind(phi)*sind(w)]

	D(5, 0) = cosd(K)*cosd(w) - sind(K)*sind(phi)*sind(w);
	D(5, 1) = sind(K)*cosd(phi)*cosd(w);
	D(5, 2) = cosd(K)*cosd(w)*sind(phi) - sind(K)*sind(w);
	//[cosd(K)*cosd(w) - sind(K)*sind(phi)*sind(w), sind(K)*cosd(phi)*cosd(w), cosd(K)*cosd(w)*sind(phi) - sind(K)*sind(w)]

	D(6, 1) = cosd(phi);
	//[0, cosd(phi),  0]

	D(7, 0) = -cosd(phi)*cosd(w);
	D(7, 1) = sind(phi)*sind(w);
	//[ -cosd(phi)*cosd(w),sind(phi)*sind(w), 0]

	D(8, 0) = -cosd(phi)*sind(w);
	D(8, 1) = -cosd(w)*sind(phi);
	//[ -cosd(phi)*sind(w),  -cosd(w)*sind(phi), 0]

	return D;
}



//Function to compute a row of A for a point
//pointNum is the number of a point that matches the point equation
MatrixXd computeAPt(int u, int numPlanes, int planeNum, int scanNum,
	double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,
	double n_xpg, double n_ypg, double n_zpg,
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj)
{
	MatrixXd A_point = MatrixXd::Zero(1, u);


	/* //Compute elements of Rbjg and RSb rotation angles

	//Rbjg
	double Rbjg_11 = cosd(K_bjg)*cosd(phi_bjg);
	double Rbjg_12_term1 = cosd(K_bjg)*sind(phi_bjg)*sind(w_bjg);
	double Rbjg_12_term2 = sind(K_bjg)*cosd(w_bjg);
	double Rbjg_13_term1 = sind(K_bjg)*sind(w_bjg);
	//double Rbjg_13_term2 = Rbjg_22_term1*sind(phi_bjg);
	double Rbjg_21 = sind(K_bjg)*cosd(phi_bjg);
	double Rbjg_22_term1 = cosd(K_bjg)*cosd(w_bjg);
	double Rbjg_22_term2 = sind(K_bjg)*sind(phi_bjg)*sind(w_bjg);
	//Rbjg_23_term1 = Rbjg_12_term2*sind(phi_bjg);
	//Rbjg_23_term2 = cosd(K_bjg)*sind(w_bjg);
	//Rbjg_31 = -sind(phi_bjg);
	double Rbjg_32 = cosd(phi_bjg)*sind(w_bjg);
	double Rbjg_33 = cosd(phi_bjg)*cosd(w_bjg);

	//RSb
	//RSb_11 = cosd(K_Sb)*cosd(phi_Sb);
	double RSb_12_term1 = cosd(K_Sb)*sind(phi_Sb)*sind(w_Sb);
	double RSb_12_term2 = sind(K_Sb)*cosd(w_Sb);
	double RSb_13_term1 = sind(K_Sb)*sind(w_Sb);
	//double RSb_13_term2 = RSb_22_term1*sind(phi_Sb);
	double RSb_21 = sind(K_Sb)*cosd(phi_Sb);
	double RSb_22_term1 = cosd(K_Sb)*cosd(w_Sb);
	double RSb_22_term2 = sind(K_Sb)*sind(phi_Sb)*sind(w_Sb);
	//RSb_23_term1 = RSb_12_term2*sind(phi_Sb);
	//RSb_23_term2 = cosd(K_Sb)*sind(w_Sb);
	//RSb_31 = -sind(phi_Sb);
	double RSb_32 = cosd(phi_Sb)*sind(w_Sb);
	double RSb_33 = cosd(phi_Sb)*cosd(w_Sb); */


	//6 Boresight parameters
	A_point(0, 0) = (n_ypg*(sind(K_bjg)*cosd(w_bjg) + cosd(K_bjg)*sind(phi_bjg)*sind(w_bjg)) + n_zpg*(sind(K_bjg)*sind(w_bjg)
		- cosd(K_bjg)*cosd(w_bjg)*sind(phi_bjg)) + n_xpg*cosd(K_bjg)*cosd(phi_bjg));


	A_point(0, 1) = (n_ypg*(cosd(K_bjg)*cosd(w_bjg) - sind(K_bjg)*sind(phi_bjg)*sind(w_bjg)) + n_zpg*(cosd(K_bjg)*sind(w_bjg)
		+ sind(K_bjg)*cosd(w_bjg)*sind(phi_bjg)) - n_xpg*sind(K_bjg)*cosd(phi_bjg));


	A_point(0, 2) = (n_xpg*sind(phi_bjg) + n_zpg*cosd(phi_bjg)*cosd(w_bjg) - n_ypg*cosd(phi_bjg)*sind(w_bjg));

	MatrixXd B_PtWrtRot = PtEqnWrtRotSb(w_bjg, phi_bjg, K_bjg, x_sj, y_sj, z_sj, n_xpg, n_ypg, n_zpg);

	MatrixXd B_RotWrtAng = RotWrtAngles(w_Sb, phi_Sb, K_Sb);

	A_point.block(0, 3, 1, 3) = B_PtWrtRot.transpose()*B_RotWrtAng;

	B_PtWrtRot.resize(0, 0);
	B_RotWrtAng.resize(0, 0);


	//4 Plane parameters
	int i = 6 + 4 * planeNum;

	A_point(0, i) = (x_bjg + x_sj*(sind(phi_bjg)*(sind(K_Sb)*sind(w_Sb) - cosd(K_Sb)*cosd(w_Sb)*sind(phi_Sb)) - sind(K_bjg)*cosd(phi_bjg)*(sind(K_Sb)*cosd(w_Sb)
		+ cosd(K_Sb)*sind(phi_Sb)*sind(w_Sb)) + cosd(K_Sb)*cosd(K_bjg)*cosd(phi_Sb)*cosd(phi_bjg)) - y_sj*(sind(K_bjg)*cosd(phi_bjg)*(cosd(K_Sb)*cosd(w_Sb)
			- sind(K_Sb)*sind(phi_Sb)*sind(w_Sb)) - sind(phi_bjg)*(cosd(K_Sb)*sind(w_Sb) + sind(K_Sb)*cosd(w_Sb)*sind(phi_Sb))
			+ cosd(K_bjg)*sind(K_Sb)*cosd(phi_Sb)*cosd(phi_bjg)) + z_sj*(cosd(K_bjg)*cosd(phi_bjg)*sind(phi_Sb) + cosd(phi_Sb)*cosd(w_Sb)*sind(phi_bjg)
				+ sind(K_bjg)*cosd(phi_Sb)*cosd(phi_bjg)*sind(w_Sb)) + z_Sjb*sind(phi_bjg) + x_Sjb*cosd(K_bjg)*cosd(phi_bjg) - y_Sjb*sind(K_bjg)*cosd(phi_bjg));

	A_point(0, i + 1) = (y_bjg + x_Sjb*(sind(K_bjg)*cosd(w_bjg) + cosd(K_bjg)*sind(phi_bjg)*sind(w_bjg)) + y_Sjb*(cosd(K_bjg)*cosd(w_bjg)
		- sind(K_bjg)*sind(phi_bjg)*sind(w_bjg)) + x_sj*((sind(K_Sb)*cosd(w_Sb) + cosd(K_Sb)*sind(phi_Sb)*sind(w_Sb))*(cosd(K_bjg)*cosd(w_bjg)
			- sind(K_bjg)*sind(phi_bjg)*sind(w_bjg)) + cosd(K_Sb)*cosd(phi_Sb)*(sind(K_bjg)*cosd(w_bjg) + cosd(K_bjg)*sind(phi_bjg)*sind(w_bjg))
			- cosd(phi_bjg)*sind(w_bjg)*(sind(K_Sb)*sind(w_Sb) - cosd(K_Sb)*cosd(w_Sb)*sind(phi_Sb))) - y_sj*(sind(K_Sb)*cosd(phi_Sb)*(sind(K_bjg)*cosd(w_bjg)
				+ cosd(K_bjg)*sind(phi_bjg)*sind(w_bjg)) - (cosd(K_Sb)*cosd(w_Sb) - sind(K_Sb)*sind(phi_Sb)*sind(w_Sb))*(cosd(K_bjg)*cosd(w_bjg)
					- sind(K_bjg)*sind(phi_bjg)*sind(w_bjg)) + cosd(phi_bjg)*sind(w_bjg)*(cosd(K_Sb)*sind(w_Sb) + sind(K_Sb)*cosd(w_Sb)*sind(phi_Sb)))
		- z_sj*(cosd(phi_Sb)*sind(w_Sb)*(cosd(K_bjg)*cosd(w_bjg) - sind(K_bjg)*sind(phi_bjg)*sind(w_bjg)) - sind(phi_Sb)*(sind(K_bjg)*cosd(w_bjg)
			+ cosd(K_bjg)*sind(phi_bjg)*sind(w_bjg)) + cosd(phi_Sb)*cosd(phi_bjg)*cosd(w_Sb)*sind(w_bjg)) - z_Sjb*cosd(phi_bjg)*sind(w_bjg));

	A_point(0, i + 2) = (z_bjg + x_Sjb*(sind(K_bjg)*sind(w_bjg) - cosd(K_bjg)*cosd(w_bjg)*sind(phi_bjg)) + y_Sjb*(cosd(K_bjg)*sind(w_bjg)
		+ sind(K_bjg)*cosd(w_bjg)*sind(phi_bjg)) + x_sj*((sind(K_Sb)*cosd(w_Sb) + cosd(K_Sb)*sind(phi_Sb)*sind(w_Sb))*(cosd(K_bjg)*sind(w_bjg)
			+ sind(K_bjg)*cosd(w_bjg)*sind(phi_bjg)) + cosd(K_Sb)*cosd(phi_Sb)*(sind(K_bjg)*sind(w_bjg) - cosd(K_bjg)*cosd(w_bjg)*sind(phi_bjg))
			+ cosd(phi_bjg)*cosd(w_bjg)*(sind(K_Sb)*sind(w_Sb) - cosd(K_Sb)*cosd(w_Sb)*sind(phi_Sb))) + y_sj*((cosd(K_Sb)*cosd(w_Sb)
				- sind(K_Sb)*sind(phi_Sb)*sind(w_Sb))*(cosd(K_bjg)*sind(w_bjg) + sind(K_bjg)*cosd(w_bjg)*sind(phi_bjg))
				- sind(K_Sb)*cosd(phi_Sb)*(sind(K_bjg)*sind(w_bjg) - cosd(K_bjg)*cosd(w_bjg)*sind(phi_bjg)) + cosd(phi_bjg)*cosd(w_bjg)*(cosd(K_Sb)*sind(w_Sb)
					+ sind(K_Sb)*cosd(w_Sb)*sind(phi_Sb))) + z_sj*(sind(phi_Sb)*(sind(K_bjg)*sind(w_bjg) - cosd(K_bjg)*cosd(w_bjg)*sind(phi_bjg))
						- cosd(phi_Sb)*sind(w_Sb)*(cosd(K_bjg)*sind(w_bjg) + sind(K_bjg)*cosd(w_bjg)*sind(phi_bjg))
						+ cosd(phi_Sb)*cosd(phi_bjg)*cosd(w_Sb)*cosd(w_bjg)) + z_Sjb*cosd(phi_bjg)*cosd(w_bjg));
	/*A_point(0, i) = 0.5;
	A_point(0, i + 1) = 0.4;
	A_point(0, i + 2) = 0.1;*/
	A_point(0, i + 3) = 1;



	//6 scan parameters
	int j = 6 + 4 * numPlanes + 6 * scanNum;

	A_point(0, j) = n_xpg;

	A_point(0, j + 1) = n_ypg;

	A_point(0, j + 2) = n_zpg;

	MatrixXd S_PtWrtRot = PtEqnWrtRotbjg(x_Sjb, y_Sjb, z_Sjb, w_Sb, phi_Sb, K_Sb, x_sj, y_sj, z_sj, n_xpg, n_ypg, n_zpg);

	MatrixXd S_RotWrtAng = RotWrtAngles(w_bjg, phi_bjg, K_bjg);

	A_point.block(0, j + 3, 1, 3) = S_PtWrtRot.transpose()*S_RotWrtAng;

	S_PtWrtRot.resize(0, 0);
	S_RotWrtAng.resize(0, 0);


	//cout << A_point << endl << endl;

	return A_point;
}



MatrixXd computewScan(double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_GPS, double y_GPS, double z_GPS,
	double w_INS, double phi_INS, double K_INS)
{
	MatrixXd w1(6, 1);

	w1(0, 0) = x_bjg - x_GPS;
	w1(1, 0) = y_bjg - y_GPS;
	w1(2, 0) = z_bjg - z_GPS;
	w1(3, 0) = w_bjg - w_INS;
	w1(4, 0) = phi_bjg - phi_INS;
	w1(5, 0) = K_bjg - K_INS;

	return w1;
}



MatrixXd computewPlane(double n_xpg, double n_ypg, double n_zpg)
{
	MatrixXd w2(1, 1);

	w2(0, 0) = pow(n_xpg, 2) + pow(n_ypg, 2) + pow(n_zpg, 2) - 1;

	return w2;
}



MatrixXd computewPt(double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,
	double n_xpg, double n_ypg, double n_zpg, double d_p,
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj)
{
	MatrixXd w3(1, 1);

	//Compute elements of Rbjg and RSb rotation angles

	//Rbjg
	/* 	double Rbjg_11 = cosd(K_bjg)*cosd(phi_bjg);
	double Rbjg_12_term1 = cosd(K_bjg)*sind(phi_bjg)*sind(w_bjg);
	double Rbjg_12_term2 = sind(K_bjg)*cosd(w_bjg);
	double Rbjg_13_term1 = sind(K_bjg)*sind(w_bjg);
	//double Rbjg_13_term2 = Rbjg_22_term1*sind(phi_bjg);
	double Rbjg_21 = sind(K_bjg)*cosd(phi_bjg);
	double Rbjg_22_term1 = cosd(K_bjg)*cosd(w_bjg);
	double Rbjg_22_term2 = sind(K_bjg)*sind(phi_bjg)*sind(w_bjg);
	//Rbjg_23_term1 = Rbjg_12_term2*sind(phi_bjg);
	//Rbjg_23_term2 = cosd(K_bjg)*sind(w_bjg);
	//Rbjg_31 = -sind(phi_bjg);
	double Rbjg_32 = cosd(phi_bjg)*sind(w_bjg);
	double Rbjg_33 = cosd(phi_bjg)*cosd(w_bjg);

	//RSb
	//RSb_11 = cosd(K_Sb)*cosd(phi_Sb);
	double RSb_12_term1 = cosd(K_Sb)*sind(phi_Sb)*sind(w_Sb);
	double RSb_12_term2 = sind(K_Sb)*cosd(w_Sb);
	double RSb_13_term1 = sind(K_Sb)*sind(w_Sb);
	//double RSb_13_term2 = RSb_22_term1*sind(phi_Sb);
	double RSb_21 = sind(K_Sb)*cosd(phi_Sb);
	double RSb_22_term1 = cosd(K_Sb)*cosd(w_Sb);
	double RSb_22_term2 = sind(K_Sb)*sind(phi_Sb)*sind(w_Sb);
	//RSb_23_term1 = RSb_12_term2*sind(phi_Sb);
	//RSb_23_term2 = cosd(K_Sb)*sind(w_Sb);
	//RSb_31 = -sind(phi_Sb);
	double RSb_32 = cosd(phi_Sb)*sind(w_Sb);
	double RSb_33 = cosd(phi_Sb)*cosd(w_Sb); */


	w3(0, 0) = (d_p + n_ypg*(y_bjg + x_Sjb*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + y_Sjb*(cos(K_bjg)*cos(w_bjg)
		- sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + x_sj*((sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*cos(w_bjg)
			- sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg))
			- cos(phi_bjg)*sin(w_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)))
		- y_sj*(sin(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - (cos(K_Sb)*cos(w_Sb)
			- sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*sin(w_Sb)
				+ sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))) - z_sj*(cos(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg))
					- sin(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*sin(w_bjg))
		- z_Sjb*cos(phi_bjg)*sin(w_bjg)) + n_xpg*(x_bjg + x_sj*(sin(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb))
			- sin(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + cos(K_Sb)*cos(K_bjg)*cos(phi_Sb)*cos(phi_bjg))
			- y_sj*(sin(K_bjg)*cos(phi_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - sin(phi_bjg)*(cos(K_Sb)*sin(w_Sb)
				+ sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(K_bjg)*sin(K_Sb)*cos(phi_Sb)*cos(phi_bjg)) + z_sj*(cos(K_bjg)*cos(phi_bjg)*sin(phi_Sb)
					+ cos(phi_Sb)*cos(w_Sb)*sin(phi_bjg) + sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*sin(w_Sb)) + z_Sjb*sin(phi_bjg)
			+ x_Sjb*cos(K_bjg)*cos(phi_bjg) - y_Sjb*sin(K_bjg)*cos(phi_bjg)) + n_zpg*(z_bjg + x_Sjb*(sin(K_bjg)*sin(w_bjg)
				- cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + y_Sjb*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))
				+ x_sj*((sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))
					+ cos(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + cos(phi_bjg)*cos(w_bjg)*(sin(K_Sb)*sin(w_Sb)
						- cos(K_Sb)*cos(w_Sb)*sin(phi_Sb))) + y_sj*((cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*sin(w_bjg)
							+ sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - sin(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg))
							+ cos(phi_bjg)*cos(w_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))) + z_sj*(sin(phi_Sb)*(sin(K_bjg)*sin(w_bjg)
								- cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - cos(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))
								+ cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*cos(w_bjg)) + z_Sjb*cos(phi_bjg)*cos(w_bjg)));


	return w3;
}





