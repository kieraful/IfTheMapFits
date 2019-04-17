#include "LS.h"

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
}



//Receives the name of a file "FileName" containing a numerical matrix, and read the matrix data into variable "m"
void Read_Mat(char *FileName, MatrixXd& m) {

	m.conservativeResize(0, 0);

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
}

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
void computeAandw(MatrixXd& A, MatrixXd& H, MatrixXd& w, MatrixXd& V, int u, int numPlanes, int numScans, int numLidPts, MatrixXd& bs_params, MatrixXd& obs_bs, MatrixXd& estscans, MatrixXd& plane_details, MatrixXd& scene_details, MatrixXd& point_details)
{
	//MatrixXd& A, MatrixXd& H, MatrixXd& w, MatrixXd& V, 

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
		A1.block(A1.rows() - 6, 0, 6, u) = computeA1(u, numPlanes, i);

		w1.conservativeResize(w1.rows() + 6, 1);
		w1.block(w1.rows() - 6, 0, 6, 1) = computew1(
			estscans(i, 0), estscans(i, 1), estscans(i, 2),
			estscans(i, 3), estscans(i, 4), estscans(i, 5),
			scene_details(i, 0), scene_details(i, 1), scene_details(i, 2),
			scene_details(i, 3), scene_details(i, 4), scene_details(i, 5));
	}

	//cout << endl << A1 << endl << endl;
	//cout << endl << w1 << endl << endl;


	for (int j = 0; j<numPlanes; j++) {
		A2.conservativeResize(A2.rows() + 1, u);
		A2.block(A2.rows() - 1, 0, 1, u) = computeA2(u, j, plane_details(j, 0), plane_details(j, 1), plane_details(j, 2));

		w2.conservativeResize(w2.rows() + 1, 1);
		w2.block(w2.rows() - 1, 0, 1, 1) = computew2(plane_details(j, 0), plane_details(j, 1), plane_details(j, 2));
	}
	//cout << A2 << endl << endl;
	//cout << endl << w2 << endl << endl;


	int planeID = 0;
	int scanID = 0;
	for (int k = 0; k<numLidPts; k++) {

		A3.conservativeResize(A3.rows() + 1, u);
		w3.conservativeResize(w3.rows() + 1, 1);

		planeID = point_details(k, 3);
		scanID = point_details(k, 4);

		A3.block(A3.rows() - 1, 0, 1, u) = computeA3(u, numPlanes, planeID, scanID,
			bs_params(0, 0), bs_params(1, 0), bs_params(2, 0),
			bs_params(3, 0), bs_params(4, 0), bs_params(5, 0),
			plane_details(planeID, 0), plane_details(planeID, 1), plane_details(planeID, 2),
			estscans(scanID, 0), estscans(scanID, 1), estscans(scanID, 2),
			estscans(scanID, 3), estscans(scanID, 4), estscans(scanID, 5),
			point_details(k, 0), point_details(k, 1), point_details(k, 2));

		w3.block(w3.rows() - 1, 0, 1, 1) = computew3(bs_params(0, 0), bs_params(1, 0), bs_params(2, 0),
			bs_params(3, 0), bs_params(4, 0), bs_params(5, 0),
			plane_details(planeID, 0), plane_details(planeID, 1), plane_details(planeID, 2), plane_details(planeID, 3),
			estscans(scanID, 0), estscans(scanID, 1), estscans(scanID, 2),
			estscans(scanID, 3), estscans(scanID, 4), estscans(scanID, 5),
			point_details(k, 0), point_details(k, 1), point_details(k, 2));
	}


	//boresight offset equations
	A.conservativeResize(A.rows() + 3, u);
	A.block(A.rows() - 3, 0, 3, u) = MatrixXd::Zero(3, u);
	A(A.rows() - 3, 0) = 1;
	A(A.rows() - 2, 1) = 1;
	A(A.rows() - 1, 2) = 1;

	w.conservativeResize(w.rows() + 3, 1);
	w.block(w.rows() - 3, 0, 3, 1) = bs_params.block(0, 0, 3, 1) - obs_bs;



	A.block(0, 0, A1.rows(), u) = A1;
	A.block(numScans * 6, 0, A3.rows(), u) = A3;

	H = A2;

	A1.resize(0, 0);
	A2.resize(0, 0);
	A3.resize(0, 0);

	w.block(0, 0, w1.rows(), 1) = w1;
	w.block(numScans * 6, 0, w3.rows(), 1) = w3;

	V = w2;


	return;
}



//Function to compute 6 rows of A for a scan
//scanNum is the number of the scan that matches the scan equations
MatrixXd computeA1(int u, int numPlanes, int scanNum)
{
	MatrixXd A1 = MatrixXd::Zero(6, u);

	int i = 6 + 4 * numPlanes + 6 * scanNum;
	int k = 0;

	for (int j = i; k<6; j++) {
		A1(k, j) = 1;
		k++;
	}

	return A1;
}



//Function to compute a row of A for a plane
//planeNum is the number of the plane that matches the plane equation
MatrixXd computeA2(int u, int planeNum, double n_xpg, double n_ypg, double n_zpg)
{
	MatrixXd A2 = MatrixXd::Zero(1, u);

	int i = 6 + 4 * planeNum;

	A2(0, i) = 2 * n_xpg;
	A2(0, i + 1) = 2 * n_ypg;
	A2(0, i + 2) = 2 * n_zpg;

	return A2;
}


//Function to compute elements of a rotation matrix, image to object rotation
MatrixXd RotMatElements(double w, double phi, double K)
{
	MatrixXd Rotw(3, 3);
	MatrixXd Rotphi(3, 3);
	MatrixXd RotK(3, 3);

	Rotw << 1, 0, 0,
		0, cos(w), sin(w),
		0, -sin(w), cos(w);

	Rotphi << cos(phi), 0, -sin(phi),
		0, 1, 0,
		sin(phi), 0, cos(phi);

	RotK << cos(K), sin(K), 0,
		-sin(K), cos(K), 0,
		0, 0, 1;

	MatrixXd RotMat = RotK*Rotphi*Rotw;
	//RotMat.transposeInPlace();

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

	D(0, 1) = -cos(K)*sin(phi);
	D(0, 2) = -sin(K)*cos(phi);
	//[0, -cos(K)*sin(phi), -sin(K)*cos(phi)]

	D(1, 0) = cos(K)*cos(w)*sin(phi) - sin(K)*sin(w);
	D(1, 1) = cos(K)*cos(phi)*sin(w);
	D(1, 2) = cos(K)*cos(w) - sin(K)*sin(phi)*sin(w);
	//[cos(K)*cos(w)*sin(phi) - sin(K)*sin(w), cos(K)*cos(phi)*sin(w), cos(K)*cos(w) - sin(K)*sin(phi)*sin(w)]

	D(2, 0) = sin(K)*cos(w) + cos(K)*sin(phi)*sin(w);
	D(2, 1) = -cos(K)*cos(phi)*cos(w);
	D(2, 2) = cos(K)*sin(w) + sin(K)*cos(w)*sin(phi);
	//[sin(K)*cos(w) + cos(K)*sin(phi)*sin(w), -cos(K)*cos(phi)*cos(w), cos(K)*sin(w) + sin(K)*cos(w)*sin(phi)]

	D(3, 1) = sin(K)*sin(phi);
	D(3, 2) = -cos(K)*cos(phi);
	//[0, sin(K)*sin(phi), -cos(K)*cos(phi)]

	D(4, 0) = -cos(K)*sin(w) - sin(K)*cos(w)*sin(phi);
	D(4, 1) = -sin(K)*cos(phi)*sin(w);
	D(4, 2) = -sin(K)*cos(w) - cos(K)*sin(phi)*sin(w);
	//[ - cos(K)*sin(w) - sin(K)*cos(w)*sin(phi), -sin(K)*cos(phi)*sin(w), - sin(K)*cos(w) - cos(K)*sin(phi)*sin(w)]

	D(5, 0) = cos(K)*cos(w) - sin(K)*sin(phi)*sin(w);
	D(5, 1) = sin(K)*cos(phi)*cos(w);
	D(5, 2) = cos(K)*cos(w)*sin(phi) - sin(K)*sin(w);
	//[cos(K)*cos(w) - sin(K)*sin(phi)*sin(w), sin(K)*cos(phi)*cos(w), cos(K)*cos(w)*sin(phi) - sin(K)*sin(w)]

	D(6, 1) = cos(phi);
	//[0, cos(phi),  0]

	D(7, 0) = -cos(phi)*cos(w);
	D(7, 1) = sin(phi)*sin(w);
	//[ -cos(phi)*cos(w),sin(phi)*sin(w), 0]

	D(8, 0) = -cos(phi)*sin(w);
	D(8, 1) = -cos(w)*sin(phi);
	//[ -cos(phi)*sin(w),  -cos(w)*sin(phi), 0]

	return D;
}



//Function to compute a row of A for a point
//pointNum is the number of a point that matches the point equation
MatrixXd computeA3(int u, int numPlanes, int planeNum, int scanNum,
	double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,
	double n_xpg, double n_ypg, double n_zpg,
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj)
{
	MatrixXd A3 = MatrixXd::Zero(1, u);


	//6 Boresight parameters
	A3(0, 0) = n_zpg*sin(phi_bjg) + n_xpg*cos(K_bjg)*cos(phi_bjg) - n_ypg*sin(K_bjg)*cos(phi_bjg);

	A3(0, 1) = n_xpg*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + n_ypg*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - n_zpg*cos(phi_bjg)*sin(w_bjg);

	A3(0, 2) = n_xpg*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + n_ypg*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + n_zpg*cos(phi_bjg)*cos(w_bjg);
	
	A3(0, 3) = -1.0*n_ypg*(y_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - sin(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))) + z_sj*(sin(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - (cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)))) - n_zpg*(y_sj*(sin(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) - cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*cos(w_bjg)) + z_sj*(cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - sin(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + cos(phi_Sb)*cos(phi_bjg)*cos(w_bjg)*sin(w_Sb))) - n_xpg*(y_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(w_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg))) - z_sj*((cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - cos(phi_Sb)*sin(w_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg))));

	A3(0, 4) = n_xpg*(x_sj*(cos(phi_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + sin(K_Sb)*sin(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - cos(K_Sb)*cos(K_bjg)*cos(phi_bjg)*sin(phi_Sb)) - z_sj*(cos(w_Sb)*sin(phi_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - sin(K_Sb)*cos(phi_Sb)*cos(w_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_Sb)*cos(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)) + y_sj*(sin(phi_Sb)*sin(w_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - sin(K_Sb)*cos(phi_Sb)*sin(w_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_Sb)*cos(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*sin(w_Sb))) - n_zpg*(z_sj*(cos(K_Sb)*cos(phi_Sb)*cos(w_Sb)*sin(phi_bjg) + cos(phi_bjg)*cos(w_Sb)*cos(w_bjg)*sin(phi_Sb) + sin(K_Sb)*cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*sin(w_bjg)) - y_sj*(cos(K_Sb)*cos(phi_Sb)*sin(phi_bjg)*sin(w_Sb) + cos(phi_bjg)*cos(w_bjg)*sin(phi_Sb)*sin(w_Sb) + sin(K_Sb)*cos(phi_Sb)*cos(phi_bjg)*sin(w_Sb)*sin(w_bjg)) + x_sj*(cos(K_Sb)*sin(phi_Sb)*sin(phi_bjg) - cos(phi_Sb)*cos(phi_bjg)*cos(w_bjg) + sin(K_Sb)*cos(phi_bjg)*sin(phi_Sb)*sin(w_bjg))) + n_ypg*(x_sj*(cos(phi_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + sin(K_Sb)*sin(phi_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_Sb)*sin(K_bjg)*cos(phi_bjg)*sin(phi_Sb)) - y_sj*(sin(K_Sb)*cos(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - sin(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + cos(K_Sb)*sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*sin(w_Sb)) + z_sj*(sin(K_Sb)*cos(phi_Sb)*cos(w_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - cos(w_Sb)*sin(phi_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + cos(K_Sb)*sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)));

	A3(0, 5) = n_zpg*(y_sj*(sin(phi_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + cos(phi_bjg)*sin(w_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb))) + z_sj*(sin(phi_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_bjg)*sin(w_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb))) - x_sj*(sin(K_Sb)*cos(phi_Sb)*sin(phi_bjg) - cos(K_Sb)*cos(phi_Sb)*cos(phi_bjg)*sin(w_bjg))) - n_ypg*(y_sj*((sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + sin(K_bjg)*cos(phi_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))) + z_sj*((sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + sin(K_bjg)*cos(phi_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))) + x_sj*(cos(K_Sb)*cos(phi_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - sin(K_Sb)*sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg))) - n_xpg*(y_sj*((sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - cos(K_bjg)*cos(phi_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))) + z_sj*((sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - cos(K_bjg)*cos(phi_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))) + x_sj*(cos(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_bjg)*sin(K_Sb)*cos(phi_Sb)*cos(phi_bjg)));


	//4 Plane parameters
	int i = 6 + 4 * planeNum;

	A3(0, i) = x_bjg + y_Sjb*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + z_Sjb*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + x_sj*(sin(phi_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - sin(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_Sb)*cos(K_bjg)*cos(phi_Sb)*cos(phi_bjg)) + z_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(w_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg))) + y_sj*((cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - cos(phi_Sb)*sin(w_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg))) + x_Sjb*cos(K_bjg)*cos(phi_bjg);

	A3(0, i + 1) = y_bjg + y_Sjb*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + z_Sjb*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - x_sj*(sin(K_Sb)*cos(phi_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - sin(phi_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + cos(K_Sb)*sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)) + z_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - sin(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))) - y_sj*(sin(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - (cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))) - x_Sjb*sin(K_bjg)*cos(phi_bjg);

	A3(0, i + 2) = z_bjg + x_sj*(cos(K_Sb)*cos(phi_Sb)*sin(phi_bjg) + cos(phi_bjg)*cos(w_bjg)*sin(phi_Sb) + sin(K_Sb)*cos(phi_Sb)*cos(phi_bjg)*sin(w_bjg)) + x_Sjb*sin(phi_bjg) + z_sj*(sin(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) - cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*cos(w_bjg)) - y_sj*(cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - sin(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + cos(phi_Sb)*cos(phi_bjg)*cos(w_bjg)*sin(w_Sb)) + z_Sjb*cos(phi_bjg)*cos(w_bjg) - y_Sjb*cos(phi_bjg)*sin(w_bjg);

	A3(0, i + 3) = 1;

	
	//6 scan parameters
	int j = 6 + 4 * numPlanes + 6 * scanNum;

	A3(0, j) = n_xpg;

	A3(0, j + 1) = n_ypg;

	A3(0, j + 2) = n_zpg;

	A3(0, j + 3) = -1.0*n_xpg*(y_Sjb*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - z_Sjb*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + z_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - cos(phi_Sb)*cos(w_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg))) + y_sj*((cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + cos(phi_Sb)*sin(w_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg))) - x_sj*(sin(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + sin(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)))) - n_ypg*(y_Sjb*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - z_Sjb*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + z_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - cos(phi_Sb)*cos(w_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg))) + y_sj*((cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + cos(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg))) - x_sj*(sin(phi_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + sin(K_Sb)*cos(phi_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)))) - n_zpg*(z_sj*(cos(phi_bjg)*cos(w_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*sin(w_bjg)) + y_sj*(cos(phi_bjg)*cos(w_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - cos(phi_Sb)*cos(phi_bjg)*sin(w_Sb)*sin(w_bjg)) + x_sj*(cos(phi_bjg)*sin(phi_Sb)*sin(w_bjg) - sin(K_Sb)*cos(phi_Sb)*cos(phi_bjg)*cos(w_bjg)) + y_Sjb*cos(phi_bjg)*cos(w_bjg) + z_Sjb*cos(phi_bjg)*sin(w_bjg));

	A3(0, j + 4) = n_ypg*(x_sj*(cos(K_Sb)*sin(K_bjg)*cos(phi_Sb)*sin(phi_bjg) + sin(K_bjg)*cos(phi_bjg)*cos(w_bjg)*sin(phi_Sb) + sin(K_Sb)*sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*sin(w_bjg)) - y_sj*(sin(K_bjg)*cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - sin(K_bjg)*sin(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*cos(w_bjg)*sin(w_Sb)) + z_sj*(sin(K_bjg)*sin(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) - sin(K_bjg)*cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*cos(w_bjg)) + x_Sjb*sin(K_bjg)*sin(phi_bjg) + z_Sjb*sin(K_bjg)*cos(phi_bjg)*cos(w_bjg) - y_Sjb*sin(K_bjg)*cos(phi_bjg)*sin(w_bjg)) - n_xpg*(x_sj*(cos(K_Sb)*cos(K_bjg)*cos(phi_Sb)*sin(phi_bjg) + cos(K_bjg)*cos(phi_bjg)*cos(w_bjg)*sin(phi_Sb) + cos(K_bjg)*sin(K_Sb)*cos(phi_Sb)*cos(phi_bjg)*sin(w_bjg)) + z_sj*(cos(K_bjg)*sin(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) - cos(K_bjg)*cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*cos(w_Sb)*cos(w_bjg)) - y_sj*(cos(K_bjg)*cos(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - cos(K_bjg)*sin(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + cos(K_bjg)*cos(phi_Sb)*cos(phi_bjg)*cos(w_bjg)*sin(w_Sb)) + x_Sjb*cos(K_bjg)*sin(phi_bjg) + z_Sjb*cos(K_bjg)*cos(phi_bjg)*cos(w_bjg) - y_Sjb*cos(K_bjg)*cos(phi_bjg)*sin(w_bjg)) + n_zpg*(x_Sjb*cos(phi_bjg) - x_sj*(cos(w_bjg)*sin(phi_Sb)*sin(phi_bjg) - cos(K_Sb)*cos(phi_Sb)*cos(phi_bjg) + sin(K_Sb)*cos(phi_Sb)*sin(phi_bjg)*sin(w_bjg)) + z_sj*(cos(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + sin(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb)) - cos(phi_Sb)*cos(w_Sb)*cos(w_bjg)*sin(phi_bjg)) + y_sj*(cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + sin(phi_bjg)*sin(w_bjg)*(cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb)) + cos(phi_Sb)*cos(w_bjg)*sin(phi_bjg)*sin(w_Sb)) - z_Sjb*cos(w_bjg)*sin(phi_bjg) + y_Sjb*sin(phi_bjg)*sin(w_bjg));

	A3(0, j + 5) = n_xpg*(y_Sjb*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + z_Sjb*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - x_sj*(sin(K_Sb)*cos(phi_Sb)*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - sin(phi_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + cos(K_Sb)*sin(K_bjg)*cos(phi_Sb)*cos(phi_bjg)) + z_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) - sin(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))) - y_sj*(sin(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - (cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(phi_Sb)*sin(w_Sb)*(cos(K_bjg)*sin(w_bjg) + sin(K_bjg)*cos(w_bjg)*sin(phi_bjg))) - x_Sjb*sin(K_bjg)*cos(phi_bjg)) - n_ypg*(y_Sjb*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + z_Sjb*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) + x_sj*(sin(phi_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg)) - sin(K_Sb)*cos(phi_Sb)*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_Sb)*cos(K_bjg)*cos(phi_Sb)*cos(phi_bjg)) + z_sj*((cos(K_Sb)*sin(w_Sb) + sin(K_Sb)*cos(w_Sb)*sin(phi_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb)) + cos(phi_Sb)*cos(w_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg))) + y_sj*((cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(phi_Sb)*sin(w_Sb))*(sin(K_bjg)*cos(w_bjg) + cos(K_bjg)*sin(phi_bjg)*sin(w_bjg)) + cos(K_bjg)*cos(phi_bjg)*(sin(K_Sb)*cos(w_Sb) + cos(K_Sb)*sin(phi_Sb)*sin(w_Sb)) - cos(phi_Sb)*sin(w_Sb)*(sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg))) + x_Sjb*cos(K_bjg)*cos(phi_bjg));


	return A3;
}



MatrixXd computew1(double x_bjg, double y_bjg, double z_bjg,
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



MatrixXd computew2(double n_xpg, double n_ypg, double n_zpg)
{
	MatrixXd w2(1, 1);

	w2(0, 0) = pow(n_xpg, 2) + pow(n_ypg, 2) + pow(n_zpg, 2) - 1;

	return w2;
}



MatrixXd computew3(double x_Sjb, double y_Sjb, double z_Sjb,
	double w_Sb, double phi_Sb, double K_Sb,
	double n_xpg, double n_ypg, double n_zpg, double d_p,
	double x_bjg, double y_bjg, double z_bjg,
	double w_bjg, double phi_bjg, double K_bjg,
	double x_sj, double y_sj, double z_sj)
{
	MatrixXd w3(1, 1);

	double cosw_Sb = cos(w_Sb);
	double cosphi_Sb = cos(phi_Sb);
	double cosK_Sb = cos(K_Sb);
	double cosw_bjg = cos(w_bjg);
	double cosphi_bjg = cos(phi_bjg);
	double cosK_bjg = cos(K_bjg);

	double sinw_Sb = sin(w_Sb);
	double sinphi_Sb = sin(phi_Sb);
	double sinK_Sb = sin(K_Sb);
	double sinw_bjg = sin(w_bjg);
	double sinphi_bjg = sin(phi_bjg);
	double sinK_bjg = sin(K_bjg);


	w3(0, 0) = (d_p + n_ypg*(y_bjg + x_Sjb*(sinK_bjg*cosw_bjg + cosK_bjg*sinphi_bjg*sinw_bjg) + y_Sjb*(cosK_bjg*cosw_bjg
		- sinK_bjg*sinphi_bjg*sinw_bjg) + x_sj*((sinK_Sb*cosw_Sb + cosK_Sb*sinphi_Sb*sinw_Sb)*(cosK_bjg*cosw_bjg
			- sinK_bjg*sinphi_bjg*sinw_bjg) + cosK_Sb*cosphi_Sb*(sinK_bjg*cosw_bjg + cosK_bjg*sinphi_bjg*sinw_bjg)
			- cosphi_bjg*sinw_bjg*(sinK_Sb*sinw_Sb - cosK_Sb*cosw_Sb*sinphi_Sb))
		- y_sj*(sinK_Sb*cosphi_Sb*(sinK_bjg*cosw_bjg + cosK_bjg*sinphi_bjg*sinw_bjg) - (cosK_Sb*cosw_Sb
			- sinK_Sb*sinphi_Sb*sinw_Sb)*(cosK_bjg*cosw_bjg - sinK_bjg*sinphi_bjg*sinw_bjg) + cosphi_bjg*sinw_bjg*(cosK_Sb*sinw_Sb
				+ sinK_Sb*cosw_Sb*sinphi_Sb)) - z_sj*(cosphi_Sb*sinw_Sb*(cosK_bjg*cosw_bjg - sinK_bjg*sinphi_bjg*sinw_bjg)
					- sinphi_Sb*(sinK_bjg*cosw_bjg + cosK_bjg*sinphi_bjg*sinw_bjg) + cosphi_Sb*cosphi_bjg*cosw_Sb*sinw_bjg)
		- z_Sjb*cosphi_bjg*sinw_bjg) + n_xpg*(x_bjg + x_sj*(sinphi_bjg*(sinK_Sb*sinw_Sb - cosK_Sb*cosw_Sb*sinphi_Sb)
			- sinK_bjg*cosphi_bjg*(sinK_Sb*cosw_Sb + cosK_Sb*sinphi_Sb*sinw_Sb) + cosK_Sb*cosK_bjg*cosphi_Sb*cosphi_bjg)
			- y_sj*(sinK_bjg*cosphi_bjg*(cosK_Sb*cosw_Sb - sinK_Sb*sinphi_Sb*sinw_Sb) - sinphi_bjg*(cosK_Sb*sinw_Sb
				+ sinK_Sb*cosw_Sb*sinphi_Sb) + cosK_bjg*sinK_Sb*cosphi_Sb*cosphi_bjg) + z_sj*(cosK_bjg*cosphi_bjg*sinphi_Sb
					+ cosphi_Sb*cosw_Sb*sinphi_bjg + sinK_bjg*cosphi_Sb*cosphi_bjg*sinw_Sb) + z_Sjb*sinphi_bjg
			+ x_Sjb*cosK_bjg*cosphi_bjg - y_Sjb*sinK_bjg*cosphi_bjg) + n_zpg*(z_bjg + x_Sjb*(sinK_bjg*sinw_bjg
				- cosK_bjg*cosw_bjg*sinphi_bjg) + y_Sjb*(cosK_bjg*sinw_bjg + sinK_bjg*cosw_bjg*sinphi_bjg)
				+ x_sj*((sinK_Sb*cosw_Sb + cosK_Sb*sinphi_Sb*sinw_Sb)*(cosK_bjg*sinw_bjg + sinK_bjg*cosw_bjg*sinphi_bjg)
					+ cosK_Sb*cosphi_Sb*(sinK_bjg*sinw_bjg - cosK_bjg*cosw_bjg*sinphi_bjg) + cosphi_bjg*cosw_bjg*(sinK_Sb*sinw_Sb
						- cosK_Sb*cosw_Sb*sinphi_Sb)) + y_sj*((cosK_Sb*cosw_Sb - sinK_Sb*sinphi_Sb*sinw_Sb)*(cosK_bjg*sinw_bjg
							+ sinK_bjg*cosw_bjg*sinphi_bjg) - sinK_Sb*cosphi_Sb*(sinK_bjg*sinw_bjg - cosK_bjg*cosw_bjg*sinphi_bjg)
							+ cosphi_bjg*cosw_bjg*(cosK_Sb*sinw_Sb + sinK_Sb*cosw_Sb*sinphi_Sb)) + z_sj*(sinphi_Sb*(sinK_bjg*sinw_bjg
								- cosK_bjg*cosw_bjg*sinphi_bjg) - cosphi_Sb*sinw_Sb*(cosK_bjg*sinw_bjg + sinK_bjg*cosw_bjg*sinphi_bjg)
								+ cosphi_Sb*cosphi_bjg*cosw_Sb*cosw_bjg) + z_Sjb*cosphi_bjg*cosw_bjg));

	//cout << endl << w3 << endl << endl;

	return w3;
}