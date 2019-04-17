#include "LS.h"



int main()
{	
	MatrixXd planeparams;
	MatrixXd scanparams; //OBSERVED
	MatrixXd pointparams;

	Read_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\planeparams.txt", planeparams);
	Read_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\scanparams.txt", scanparams);
	Read_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\pointparams.txt", pointparams);

	//Read_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\Plane_Details.txt", planeparams);
	//Read_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\Scene_Details.txt", scanparams);
	//Read_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\Point_Details.txt", pointparams);

	//Convert vector of vectors to MatrixXd for plane, scene, and lidar point parameters
	/*vec2mat(plane_details, planeparams, 4);
	vec2mat(scene_details, scanparams, 6);
	vec2mat(point_details, pointparams, 5);*/

	//Write_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\planeparams.txt", planeparams, 5);
	//Write_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\scanparams.txt", scanparams, 5);
	//Write_Mat("C:\\Users\\Esther\\Desktop\\500\\data\\pointparams.txt", pointparams, 5);

	//Get number of planes, scans, lidar points
	int numPlanes = planeparams.rows();
	int numScans = scanparams.rows();
	int numLidPts = pointparams.rows();


	for (int i = 0; i < numScans; i++)
	{
		scanparams(i, 3) *= PI / 180.0;
		scanparams(i, 4) *= PI / 180.0;
		scanparams(i, 5) *= PI / 180.0;
	}

	//total number of unknowns
	int u = 6 + 4 * numPlanes + 6 * numScans;
	int n = numLidPts*3 + numScans*6;

	//Initial boresight parameters set to zero
	//MatrixXd bs_params = MatrixXd::Zero(6, 1);

	//MatrixXd bs_params = MatrixXd::Ones(6, 1);
	//bs_params = bs_params*0.1;

	//Assuming rotations are image to object
	MatrixXd bs_params(6, 1);
	bs_params(0, 0) = -0.661;
	bs_params(1, 0) = -0.272;
	bs_params(2, 0) = -1.4158;
	bs_params(3, 0) = PI;
	bs_params(4, 0) = 0.0;
	bs_params(5, 0) = -1.0*PI/2;

	MatrixXd obs_bs = MatrixXd::Zero(3, 1);
	//MatrixXd obs_bs = MatrixXd::Ones(3, 1);
	//obs_bs = obs_bs*0.1;

	//obs_bs(0, 0) = -0.661;
	//obs_bs(1, 0) = -0.272;
	//obs_bs(2, 0) = -1.4158;

	//Initial estimated scan EOPs set to zero
	MatrixXd estscans = scanparams;
	
	//Initialize matrices for LS
	MatrixXd A(numScans * 6 + numLidPts, u);
	MatrixXd A_full(numScans * 6 + numPlanes + numLidPts, u);
	MatrixXd N(0, 0);
	MatrixXd w(numScans * 6 + numLidPts, 1);
	MatrixXd w_full(numScans * 6 + numPlanes + numLidPts, 1);
	MatrixXd U;
	MatrixXd H(numPlanes, u);
	MatrixXd V(numPlanes, 1);
	MatrixXd B(u + numPlanes, u + numPlanes);
	//MatrixXd B(0,0);
	MatrixXd C(u + numPlanes, 1);
	MatrixXd Y = MatrixXd::Ones(u + numPlanes, 1);
	MatrixXd v;
	MatrixXd apostvf;
	MatrixXd Cxhat;

	//MatrixXd M;
	//MatrixXd Minv;
	//MatrixXd delta;
	//MatrixXd k;

	//Change hardcoded values to actual values******************************************************
	//Set P. Clo is a diagonal matrix so inverse of Clo is a diagonal matrix with each diagonal element inversed.
	//MatrixXd P_full = MatrixXd::Identity(numScans * 6 + numPlanes + numLidPts, numScans * 6 + numPlanes + numLidPts);
	//for (int i = 0; i < (numScans * 6); i = i + 6)
	//{
	//	P_full(i, i) = 1 / pow(0.003, 2);
	//	P_full(i + 1, i + 1) = 1 / pow(0.003, 2);
	//	P_full(i + 2, i + 2) = 1 / pow(0.004, 2);
	//	P_full(i + 3, i + 3) = 1 / pow(0.00069 * 180 / PI, 2);
	//	P_full(i + 4, i + 4) = 1 / pow(0.00070 * 180 / PI, 2);
	//	P_full(i + 5, i + 5) = 1 / pow((0.00466 + 0.00289 + 0.00339) / 3 * 180 / PI, 2);
	//}
	//P_full.block(numScans * 6, numScans * 6, numPlanes, numPlanes) = MatrixXd::Identity(numPlanes, numPlanes) * 1 / sqrt(pow(0.003, 2) + pow(0.003, 2) + pow(0.004, 2));
	//P_full.block(numScans * 6 + numPlanes, numScans * 6 + numPlanes, numLidPts, numLidPts) = MatrixXd::Identity(numLidPts, numLidPts) * 1 / 0.03;

	//MatrixXd P(numScans * 6 + numLidPts, numScans * 6 + numLidPts);
	//P.block(0, 0, numScans * 6, numScans * 6) = P_full.block(0, 0, numScans * 6, numScans * 6);
	//P.block(numScans * 6, numScans * 6, numLidPts, numLidPts) = P_full.block(numScans * 6 + numPlanes, numScans * 6 + numPlanes, numLidPts, numLidPts);

	MatrixXd P = MatrixXd::Identity(numScans * 6 + numLidPts, numScans * 6 + numLidPts); // *0.1;
	//MatrixXd Cl_full = MatrixXd::Identity(n,n); // *0.1;
	//P = P*pow(0.1, 2);


	int iter = 1;
	//double mean_Y = 1;
	//mean_Y > 0.000001
	//cout << "Y: " << Y.block(0, 0, 6, 1);

	//(abs(Y.maxCoeff()) > 0.001 || abs(Y.minCoeff()) > 0.001)
	while(abs(Y.block(0,0,6,1).mean())> 0.00001)
	{
		computeAandw(A, H, w, V, u, numPlanes, numScans, numLidPts, bs_params, obs_bs, estscans, planeparams, scanparams, pointparams);

		Write_Mat("A.txt", A, 15);
		Write_Mat("H.txt", H, 15);
		Write_Mat("w.txt", w, 3);
		Write_Mat("V.txt", V, 6);

		N = (A.transpose())* A;
		Write_Mat("N.txt", N, 15);

		U = -1.0*(A.transpose())* w;

		B.block(0, 0, u, u) = N;
		B.block(u, 0, numPlanes, u) = H;
		B.block(0, u, u, numPlanes) = (H.transpose());
		B.block(u, u, numPlanes, numPlanes) = MatrixXd::Zero(numPlanes, numPlanes);
		Write_Mat("B.txt", B, 15);

		C.block(0, 0, u, 1) = U;
		C.block(u, 0, numPlanes, 1) = -V;
		Write_Mat("C.txt", C, 6);

		//FullPivLU<MatrixXd> B_decomp(B);
		//cout << "Is B invertible? " << B_decomp.isInvertible() << endl;
		//MatrixXd Binv=B.inverse();
		//Write_Mat("Binv.txt", Binv, 6);

		Y = B.ldlt().solve( C );
		Write_Mat("Y.txt", Y, 6);


		//Update unknowns

		bs_params = bs_params + Y.block(0, 0, 6, 1);

		int Yindex = 6;
		for (int i = 0; i < planeparams.rows(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				planeparams(i, j) += Y(Yindex, 0);
				Yindex++;
			}
		}

		for (int i = 0; i < estscans.rows(); i++)
		{
			for (int j = 0; j < 6; j++)
			{
				estscans(i, j) += Y(Yindex, 0);
				Yindex++;
			}
		}
		
		
		Write_Mat("boresight.txt", bs_params, 3);
		Write_Mat("estscanEOPs.txt", estscans, 3);	
		cout << "Iteration " << iter << endl;
		cout << Y.block(0, 0, 6, 1).mean() << endl;
		iter++;
	}

	//Calculate obs residuals, aposteriori variance factor and Cxhat
	v = A * Y.block(0, 0, u, 1) + w;
	apostvf = v.transpose() * v / (n - u);

	FullPivLU<MatrixXd> lu_decomp(N);
	cout << "Is N invertible? " << lu_decomp.isInvertible() << endl;
	Cxhat = apostvf(0, 0) * N.inverse();
	Write_Mat("Cxhat.txt", Cxhat, 3);

	/*P = (A_full*Cxhat*A_full.transposeInPlace());

	FullPivLU<MatrixXd> P_decomp(P);
	cout << "Is Clhat invertible? " << lu_decomp.isInvertible() << endl;

	P = (A_full*Cxhat*A_full.transposeInPlace()).inverse();
	Write_Mat("Pupdate.txt", P, 3);*/
	

	return 0;
}