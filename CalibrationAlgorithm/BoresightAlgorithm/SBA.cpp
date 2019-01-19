/*/////////////////////////////////////////////////////////////////
Copyright (c) 2016, Mozhdeh Shahbazi
All rights reserved.

Additionally, these source codes are legally licensed under
a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. 
This is not a Free Culture License.

* Adaptations of this work cannot be shared
* Commercial uses of this work are not allowed
* Using the source code must retain this license and copy-right notice

For more details see:
creativecommons.org/licenses/by-nc-nd/4.0/

Disclaimer of Warranties and Limitation of Liability:

Unless otherwise separately undertaken by the Licensor, to the extent possible, 
the Licensor offers the Licensed Material as-is and as-available, 
and makes no representations or warranties of any kind concerning the Licensed Material, 
whether express, implied, statutory, or other. This includes, without limitation, warranties of title, 
merchantability, fitness for a particular purpose, non-infringement, absence of latent or other defects, 
accuracy, or the presence or absence of errors, whether or not known or discoverable. 
Where disclaimers of warranties are not allowed in full or in part, this disclaimer may not apply to You.
To the extent possible, in no event will the Licensor be liable to You on any legal theory 
(including, without limitation, negligence) or otherwise for any direct, special, indirect, incidental, 
consequential, punitive, exemplary, or other losses, costs, expenses, or damages arising out of 
this Public License or use of the Licensed Material, even if the Licensor has been advised of the possibility 
of such losses, costs, expenses, or damages. Where a limitation of liability is not allowed in full or in part, 
this limitation may not apply to You.
The disclaimer of warranties and limitation of liability provided above shall be interpreted 
in a manner that, to the extent possible, most closely approximates an absolute disclaimer 
and waiver of all liability.
////////////////////////////////////////////////////////////////// */


#include "SBA_NoCal.h"


/*/////////////////////////////////////////////////////////////////

EDITED BY: IFTHEMAPFITS ENGO500 TEAM for Boresight Calibration. With express permission from Mozhdeh Shahbazi

////////////////////////////////////////////////////////////////// */

void Perform_BA(VectorXi & Point_Labels,
	Matrixdby2 & Obspixels_p, MatrixXd & Sigmaxy_p, Matrixdby2i & Index_p_i,
	VectorXi & Adress_i_p, VectorXi & Adress_p_i,
	VectorXi & Index_i_p, CameraParam& camera_params, MatrixXd& ground, MatrixXd delta_trans, double IO_scalfact, MatrixXd& Sigma_Xcap_C, MatrixXd& Sigma_Xcap_G) {


	//some parameters to set
	int Maxiter = 16; //one of the conditions to stop BA
	double ErrTol = 1e-15;
	double maxcorrection = 1e20;
	double Sigma_02_cap = 1.0;
	


	int n_points=Adress_p_i.size()-1; //number of tie points
	
	
	MatrixXi FrequencyP;
	FrequencyP.setZero(n_points, 1);
	for (int i = 0; i<n_points; i++) {
		FrequencyP(i, 0) = (Adress_p_i(i + 1) - Adress_p_i(i));
	}
	

	cout << "Average FrequencyP = " << FrequencyP.array().mean() << endl;
	cout << "Min FrequencyP = " << FrequencyP.array().minCoeff() << endl;


	int n_cams = Adress_i_p.size() - 1;//number of images in this network

	
	MatrixXi DensityI;
	DensityI.setZero(n_cams, 1);
	for (int j = 0; j<n_cams; j++) {
		DensityI(j, 0) = (Adress_i_p(j + 1) - Adress_i_p(j));
	}

	cout << "Average DensityI = " << DensityI.array().mean() << endl;
	cout << "Min DensityI = " << DensityI.array().minCoeff() << endl;


	//estimations of 3D coordinates of object points
	ground.setZero(n_points, 3); //X Y Z


	Matrixdb3 Bundestim_ori;
	Bundestim_ori.resize(n_cams * 5, 3);
	Bundestim_ori = camera_params.Bundestim_ori;
	for (int i = 0; i<n_cams; i++) {
		Bundestim_ori.row(i * 5 + 3) = Bundestim_ori.row(i * 5 + 3) - delta_trans;// Xo Yo Zo
	}

	Matrix3b3 ACalmat, KCalmat;
	Matrix<double, 3, 1> tmp31;
	
	//Approximate values of ground coordinates
	for (int i = 0; i<n_points; i++) {
		int eqr = 0;

		MatrixXd A_mat;
		MatrixXd L_mat;
		Matrix3b3 r;
		Matrix3b4 RT;
		Matrix3b4 Pmat1;

		
		A_mat.setZero(FrequencyP(i, 0) * 2, 3);
		L_mat.setZero(FrequencyP(i, 0) * 2, 1);
		
		int strt = Adress_p_i(i);
		int nd = Adress_p_i(i + 1) - 1;

		for (int p = strt; p <= nd; p++) {

			int im = Index_p_i(p, 1);
			int sens_im = (int)camera_params.sensor_ID(im, 0); //in your case, this is always 0
			int model_im = (int)camera_params.model_ID(im, 0); //in your case, this is always 0
		
			ACalmat << camera_params.PS_0(sens_im, 0), 0, (.5 - camera_params.Cn_0(sens_im, 0) / 2)*camera_params.PS_0(sens_im, 0) - camera_params.xpp_0(sens_im, 0),
					0, -camera_params.PS_0(sens_im, 0), -(.5 - camera_params.Rn_0(sens_im, 0) / 2)*camera_params.PS_0(sens_im, 0) - camera_params.ypp_0(sens_im, 0),
					0, 0, -camera_params.f_l_0(sens_im, 0);
			KCalmat = (ACalmat.inverse());
			KCalmat = KCalmat / KCalmat(2, 2); //Normalizing camera calibration matrix
			


			r = Bundestim_ori.block(im * 5 + 0, 0, 3, 3);
			RT.block(0, 0, 3, 3) = r;
			RT.block(0, 3, 3, 1) = -r * ((Bundestim_ori.row(im * 5 + 3)).transpose());
			Pmat1 = KCalmat * RT; //Camera perspective projection matrix

			//Image observations in sensor coordinate system
			double x1 = Obspixels_p(p, 0);
			double y1 = Obspixels_p(p, 1); 

			
			A_mat.row(eqr) << Pmat1(0, 0) - x1 * Pmat1(2, 0), Pmat1(0, 1) - x1 * Pmat1(2, 1), Pmat1(0, 2) - x1 * Pmat1(2, 2);

			A_mat.row(eqr + 1) << Pmat1(1, 0) - y1 * Pmat1(2, 0), Pmat1(1, 1) - y1 * Pmat1(2, 1), Pmat1(1, 2) - y1 * Pmat1(2, 2);

			L_mat.row(eqr) << -Pmat1(0, 3) + Pmat1(2, 3)*x1;
			L_mat.row(eqr + 1) << -Pmat1(1, 3) + Pmat1(2, 3)*y1;

			eqr = eqr + 2;
			

		}


		tmp31 = (A_mat.transpose()*A_mat).lu().solve(A_mat.transpose()*L_mat);
		ground.row(i) << tmp31(0, 0), tmp31(1, 0), tmp31(2, 0);
	}


	Write_Mat("..txt",ground,3);

	//This scale factor can help becoming ill-configured in design matrix
	MatrixXd PS_0 = IO_scalfact * camera_params.PS_0;
	MatrixXd f_l_0 = IO_scalfact * camera_params.f_l_0;
	MatrixXd xpp_0 = IO_scalfact * camera_params.xpp_0;
	MatrixXd ypp_0 = IO_scalfact * camera_params.ypp_0;
	MatrixXd Cn_0 = camera_params.Cn_0;
	MatrixXd Rn_0 = camera_params.Rn_0;
	MatrixXd K1_0 = (1 / pow(IO_scalfact, (int)2))*camera_params.K1_0;
	MatrixXd K2_0 = (1 / pow(IO_scalfact, (int)4))*camera_params.K2_0;
	MatrixXd K3_0 = (1 / pow(IO_scalfact, (int)6))*camera_params.K3_0;
	MatrixXd P1_0 = (1 / IO_scalfact)*camera_params.P1_0;
	MatrixXd P2_0 = (1 / IO_scalfact)*camera_params.P2_0;
	MatrixXd S1_0 = camera_params.S1_0;
	MatrixXd S2_0 = camera_params.S2_0;
	MatrixBools Calibrated_0 = camera_params.Calibrated_0;



	int n_eq=Obspixels_p.rows(); //number of observed image points (2 times of n_eq is the total number of collinearity equations we have in the system)
	int u_o=n_cams*6; //number of image orientations unknowns: Omega Fi Kapa Xo Yo Zo
	int u_g=n_points*3; // number of ground coordinates unknowns: X Y Z
	
    //Unknowns
	int ioparams = 10;
	int u_M = ioparams; 
	int eq_u = u_M + n_cams * 6 + n_points * 3;

		
	int Loopcounter=0;

	while (Loopcounter <= Maxiter && maxcorrection > ErrTol) {
		Loopcounter = Loopcounter + 1;

		cout << endl << "Loopcounter=" << Loopcounter << endl;

		//Dij=augmented collinearity equations
		MatrixXd AM_ij;
		AM_ij.setZero(2, u_M*n_eq); //d(Dij)/d(M)
		

		MatrixXd AC_ij;
		AC_ij.setZero(2, 6 * n_eq); //d(Dij)/d(Cj) 
		MatrixXd AG_ij;
		AG_ij.setZero(2, 3 * n_eq); //d(Dij)/d(Gi)


		MatrixXd B_ij;
		B_ij.setZero(2, 2 * n_eq); //d(Dij)/d(xyij)
		MatrixXd D_ij;
		D_ij.setZero(2, 2 * n_eq); //inv(B_ij*Pvinv_ij*B_ij')
		MatrixXd Zy_ij;
		Zy_ij.setZero(2, 2 * n_eq); //Pvinv_ij*B_ij'
		MatrixXd W1_ij;
		W1_ij.setZero(2, 1 * n_eq); //Dij evaluated at X0

		int eq_r = 0;

		int dddd = 0;
		for (int i = 0; i<n_points; i++) { //point number
			
			int strt = Adress_p_i(i);
			int nd = Adress_p_i(i + 1) - 1;

			for (int p = strt; p <= nd; p++) { //p=ij=position of this particular observation (point i observed on image j) at matrix of observations Obspixels_p

				int j = Index_p_i(p, 1); //image number

				int sens_j = (int)camera_params.sensor_ID(j, 0); //sensor index; in your case always 0 

				int model_j = (int)camera_params.model_ID(j, 0); //model index; in your case always 0

				double Oemga0 = Bundestim_ori(j * 5 + 4, 0);
				double Phi0 = Bundestim_ori(j * 5 + 4, 1);
				double Kappa0 = Bundestim_ori(j * 5 + 4, 2);

				double Xo0 = Bundestim_ori(j * 5 + 3, 0);
				double Yo0 = Bundestim_ori(j * 5 + 3, 1);
				double Zo0 = Bundestim_ori(j * 5 + 3, 2);

				Matrix3b3 mat0 = Bundestim_ori.block(j * 5, 0, 3, 3);
				Matrix<double, 1, 3> M01 = mat0.row(0);
				Matrix<double, 1, 3> M02 = mat0.row(1);
				Matrix<double, 1, 3> M03 = mat0.row(2);

				Matrix<double, 1, 3> Xbar;
				Xbar << ground(i, 0) - Xo0, ground(i, 1) - Yo0, ground(i, 2) - Zo0;
				double q = M03(0, 0)*Xbar(0, 0) + M03(0, 1)*Xbar(0, 1) + M03(0, 2)*Xbar(0, 2); //in lectures we showd this by U3_ij
				double q2 = q * q;
				double s = M02(0, 0)*Xbar(0, 0) + M02(0, 1)*Xbar(0, 1) + M02(0, 2)*Xbar(0, 2);//in lectures we showd this by U2_ij	
				double r = M01(0, 0)*Xbar(0, 0) + M01(0, 1)*Xbar(0, 1) + M01(0, 2)*Xbar(0, 2);//in lectures we showd this by U1_ij

				double xd = Obspixels_p(p, 0);
				double yd = Obspixels_p(p, 1);

				Matrix<double, 1, 2> bb1;
				Matrix<double, 1, 6> ao1;
				Matrix<double, 1, 3> ag1;
				MatrixXd am1;
				Matrix<double, 1, 6> ao2;
				Matrix<double, 1, 3> ag2;
				Matrix<double, 1, 2> bb2;
				MatrixXd am2;
				double w1, w2;

				if (model_j == 0) {
					xd = (Obspixels_p(p, 0) + 0.5 - Cn_0(sens_j, 0) / 2)*PS_0(sens_j, 0);
					yd = -(Obspixels_p(p, 1) + 0.5 - Rn_0(sens_j, 0) / 2)*PS_0(sens_j, 0);
					double xdd = (xd - xpp_0(sens_j, 0));
					double ydd = (yd - ypp_0(sens_j, 0));
					double r2 = xdd * xdd + ydd * ydd;

					double r2prime_xpp = -2 * (xd - xpp_0(sens_j, 0));
					double r2prime_ypp = -2 * (yd - ypp_0(sens_j, 0));

					double b11 = (f_l_0(sens_j, 0) / (q2))*(r*(-mat0(2, 2)*Xbar(0, 1) + mat0(2, 1)*Xbar(0, 2)) - q * (-mat0(0, 2)*Xbar(0, 1) + mat0(0, 1)*Xbar(0, 2)));
					double b12 = (f_l_0(sens_j, 0) / (q2))*(r*(cos(Phi0)*Xbar(0, 0) + sin(Oemga0)*sin(Phi0)*Xbar(0, 1) - cos(Oemga0)*sin(Phi0)*Xbar(0, 2)) - q * (-sin(Phi0)*cos(Kappa0)*Xbar(0, 0) + sin(Oemga0)*cos(Phi0)*cos(Kappa0)*Xbar(0, 1) - cos(Oemga0)*cos(Phi0)*cos(Kappa0)*Xbar(0, 2)));
					double b13 = (-f_l_0(sens_j, 0) / (q))*(mat0(1, 0)*Xbar(0, 0) + mat0(1, 1)*Xbar(0, 1) + mat0(1, 2)*Xbar(0, 2));
					double b14 = (f_l_0(sens_j, 0) / (q2))*(r*mat0(2, 0) - q * mat0(0, 0));
					double b15 = (f_l_0(sens_j, 0) / (q2))*(r*mat0(2, 1) - q * mat0(0, 1));
					double b16 = (f_l_0(sens_j, 0) / (q2))*(r*mat0(2, 2) - q * mat0(0, 2));

					w1 = (xdd + xdd * (K1_0(sens_j, 0)*r2 + K2_0(sens_j, 0) * r2*r2 + K3_0(sens_j, 0) * r2*r2*r2) + P1_0(sens_j, 0) * (r2 + 2 * xdd*xdd) + 2 * P2_0(sens_j, 0)*xdd*ydd + S1_0(sens_j, 0) * xdd + S2_0(sens_j, 0) * ydd) + f_l_0(sens_j, 0) * r / q; //F0_ij

					ao1.setZero(1, 6);
					ao1 << -b11, -b12, -b13, b14, b15, b16; //derivatives of the first collineairty equation w.r.t Omega_j, Phi_j, Kappa_j, Xo_j, Yo_j, Zo_j

					ag1.setZero(1, 3);
					ag1 << -b14, -b15, -b16;//w.r.t Xi Yi Zi


					double rondK1 = xdd * r2;  //w.r.t K1
					double rondK2 = xdd * r2*r2;  //w.r.t K2
					double rondK3 = xdd * r2*r2*r2;  //w.r.t K3
					double rondP1 = r2 + 2 * xdd*xdd;  //w.r.t P1
					double rondP2 = 2 * xdd*ydd;  //w.r.t P2
					double rondS1 = xdd;  //derivation of F_xij w.r.t S1
					double rondS2 = ydd;  //derivation of F_xij w.r.t S2

					double rondcxpp = -1 - 1 * (K1_0(sens_j, 0)*r2 + K2_0(sens_j, 0) * (r2*r2) + K3_0(sens_j, 0) * (r2*r2*r2)) + (xd - xpp_0(sens_j, 0)) * 1 * (K1_0(sens_j, 0)*r2prime_xpp + K2_0(sens_j, 0) * 2 * r2*r2prime_xpp + K3_0(sens_j, 0) * 3 * ((r2*r2))*r2prime_xpp)
						+ P1_0(sens_j, 0) * r2prime_xpp - 2 * P1_0(sens_j, 0) * 2 * (xd - xpp_0(sens_j, 0)) - 2 * P2_0(sens_j, 0) * 1 * ydd -S1_0(sens_j, 0); //w.r.t xpp; %w.r.t cxpp

					double rondcypp = (xd - xpp_0(sens_j, 0)) * 1 * (K1_0(sens_j, 0)*r2prime_ypp + K2_0(sens_j, 0) * 2 * r2*r2prime_ypp + K3_0(sens_j, 0) * 3 * ((r2*r2))*r2prime_ypp)
						+ P1_0(sens_j, 0) * r2prime_ypp - 2 * P2_0(sens_j, 0)*xdd * 1-S2_0(sens_j, 0); //w.r.t cypp

					double rondf = r / q; //w.r.t f

					double rondx = -rondcxpp * PS_0(sens_j, 0); //w.r.t obs_pixX in sensor coordinate system,pixels
					double rondy = -rondcypp * PS_0(sens_j, 0); //w.r.t obs_pixY in sensor coordinate system,pixels


					am1.setZero(1, u_M);
					am1<< rondK1, rondK2, rondK3, rondP1, rondP2, rondS1, rondS2, rondcxpp, rondcypp, rondf;


					bb1.setZero(1, 2);
					bb1 << rondx, rondy;


					double b21 = (f_l_0(sens_j, 0) / (q2))*(s*(-mat0(2, 2)*Xbar(0, 1) + mat0(2, 1)*Xbar(0, 2)) - q * (-mat0(1, 2)*Xbar(0, 1) + mat0(1, 1)*Xbar(0, 2)));
					double b22 = (f_l_0(sens_j, 0) / (q2))*(s*(cos(Phi0)*Xbar(0, 0) + sin(Oemga0)*sin(Phi0)*Xbar(0, 1) - cos(Oemga0)*sin(Phi0)*Xbar(0, 2)) - q * (sin(Phi0)*sin(Kappa0)*Xbar(0, 0) - sin(Oemga0)*cos(Phi0)*sin(Kappa0)*Xbar(0, 1) + cos(Oemga0)*cos(Phi0)*sin(Kappa0)*Xbar(0, 2)));
					double b23 = (f_l_0(sens_j, 0) / (q))*(mat0(0, 0)*Xbar(0, 0) + mat0(0, 1)*Xbar(0, 1) + mat0(0, 2)*Xbar(0, 2));
					double b24 = (f_l_0(sens_j, 0) / (q2))*(s*mat0(2, 0) - q * mat0(1, 0));
					double b25 = (f_l_0(sens_j, 0) / (q2))*(s*mat0(2, 1) - q * mat0(1, 1));
					double b26 = (f_l_0(sens_j, 0) / (q2))*(s*mat0(2, 2) - q * mat0(1, 2));

					w2 = (ydd + ydd * (K1_0(sens_j, 0)*r2 + K2_0(sens_j, 0) * (r2*r2) + K3_0(sens_j, 0) * (r2*r2*r2)) + P2_0(sens_j, 0) * (r2 + 2 * ydd*ydd) + 2 * P1_0(sens_j, 0)*xdd*ydd) + f_l_0(sens_j, 0) * s / q; //F_yij(X0) 

					ao2.setZero(1, 6);
					ao2 << -b21, -b22, -b23, b24, b25, b26; // derivatives of the second collineairty equation w.r.t Omega Phi Kappa Xc Yc Zc_j

					ag2.setZero(1, 3);
					ag2 << -b24, -b25, -b26;//w.r.t Xi Yi Zi

					rondK1 = ydd * r2;  //w.r.t K1
					rondK2 = ydd * (r2*r2);  //w.r.t K2
					rondK3 = ydd * (r2*r2*r2);  //w.r.t K3
					rondP1 = 2 * xdd*ydd;  //w.r.t P1
					rondP2 = r2 + 2 * ydd*ydd;  //w.r.t P2
					rondS1 = 0;  //derivation of F_yij w.r.t S1
					rondS2 = 0;  //derivation of F_yij w.r.t S2

					rondcxpp = (yd - ypp_0(sens_j, 0)) * 1 * (K1_0(sens_j, 0)*r2prime_xpp + K2_0(sens_j, 0) * 2 * r2*r2prime_xpp + K3_0(sens_j, 0) * 3 * ((r2*r2))*r2prime_xpp)
						+ P2_0(sens_j, 0) * r2prime_xpp - 2 * P1_0(sens_j, 0)*ydd * 1; //w.r.t cxpp

					rondcypp = -1 - 1 * (K1_0(sens_j, 0)*r2 + K2_0(sens_j, 0) * (r2*r2) + K3_0(sens_j, 0) * (r2*r2*r2)) + (yd - ypp_0(sens_j, 0)) * 1 * (K1_0(sens_j, 0)*r2prime_ypp + K2_0(sens_j, 0) * 2 * r2*r2prime_ypp + K3_0(sens_j, 0) * 3 * ((r2*r2))*r2prime_ypp)
						+ P2_0(sens_j, 0) * r2prime_ypp - 2 * P2_0(sens_j, 0) * 2 * (yd - ypp_0(sens_j, 0)) - 2 * P1_0(sens_j, 0) * 1 * xdd; //w.r.t cypp

					rondf = s / q; //w.r.t f;

					rondx = -rondcxpp * PS_0(sens_j, 0); //w.r.t obs_pixX in sensor coordinate system,pixels
					rondy = -rondcypp * PS_0(sens_j, 0); //w.r.t obs_pixY in sensor coordinate system,pixels


					bb2.setZero(1, 2);
					bb2 << rondx, rondy;


					am2.setZero(1, u_M);
					am2 << rondK1, rondK2, rondK3, rondP1, rondP2, rondS1, rondS2, rondcxpp, rondcypp, rondf;

					
				}

				W1_ij(0, p) = w1; W1_ij(1, p) = w2;

				AC_ij.block(0, 6 * p, 1, 6) = ao1; AC_ij.block(1, 6 * p, 1, 6) = ao2;

				AG_ij.block(0, 3 * p, 1, 3) = ag1; AG_ij.block(1, 3 * p, 1, 3) = ag2;

				AM_ij.block(0, u_M*p, 1, u_M) = am1;
				AM_ij.block(1, u_M*p, 1, u_M) = am2;
				
				B_ij.block(0, 2 * p, 1, 2) = bb1; B_ij.block(1, 2 * p, 1, 2) = bb2;

				Matrix<double, 2, 2> temp = ((Sigmaxy_p(p, 0)*Sigmaxy_p(p, 0))*MatrixXd::Identity(2, 2))*(B_ij.block(0, 2 * p, 2, 2).transpose());
				Zy_ij.block(0, 2 * p, 2, 2) = temp;
				temp = (B_ij.block(0, 2 * p, 2, 2)*Zy_ij.block(0, 2 * p, 2, 2));

				//instead of using .inverse() from Eigen
				Matrix<double, 2, 2> temp2;
				double detr = 1 / (temp(0, 0)*temp(1, 1) - temp(1, 0)*temp(0, 1));
				temp2 << detr * temp(1, 1), -detr * temp(0, 1),
					-detr * temp(1, 0), detr*temp(0, 0);
				D_ij.block(0, 2 * p, 2, 2) = temp2;

			}
		}

		
		MatrixXd NMM;
		if (u_M > 0) {
			NMM.setZero(u_M, u_M);  //sigma_over_all(AM_ij'*D_ij*AM_ij)
		}
		MatrixXd UM;
		if (u_M > 0) {
			UM.setZero(u_M, 1); //sigma_over_all(AM_ij'*D_ij*W1_ij)
		}

		MatrixXd tmp;
		if (u_M > 0) {
			for (int p = 0; p < n_eq; p++) {
				int j = Index_p_i(p, 1); //image number
				tmp = (AM_ij.block(0, p*u_M, 2, u_M).transpose())*D_ij.block(0, 2 * p, 2, 2);
				NMM = NMM + tmp * AM_ij.block(0, p*u_M, 2, u_M);
				UM = UM + (tmp*W1_ij.block(0, p * 1, 2, 1));
			}
			UM = -UM;
		}


		MatrixXd NMC_j;
		if (u_M > 0) {
			NMC_j.setZero(u_M, 6 * n_cams);//sigma_over_i(AM_ij'*D_ij*AC_ij)
		}
		MatrixXd NMCt_j;
		if (u_M > 0) {
			NMCt_j.setZero(6, u_M*n_cams);//sigma_over_i(AC_ij'*D_ij*AM_ij)
		}


		MatrixXd UC_j;
		UC_j.setZero(6, 1 * n_cams);//sum_over_i(AC_ij'*P_ij*W1_ij)

		MatrixXd NCC_j;
		NCC_j.setZero(6, 6 * n_cams);//sum_over_i(AC_ij'*P_ij*AC_ij)

		for (int j = 0; j<n_cams; j++) {
			int strt = Adress_i_p(j);
			int nd = Adress_i_p(j + 1) - 1;

			for (int indx = strt; indx <= nd; indx++) { //sum over i
				int p = Index_i_p(indx);
				if (u_M > 0) {
					NMC_j.block(0, j * 6, u_M, 6) = NMC_j.block(0, j * 6, u_M, 6) + (AM_ij.block(0, p*u_M, 2, u_M).transpose())*D_ij.block(0, 2 * p, 2, 2)*AC_ij.block(0, 6 * p, 2, 6);
					NMCt_j.block(0, u_M*j, 6, u_M) = NMCt_j.block(0, u_M*j, 6, u_M) + (AC_ij.block(0, 6 * p, 2, 6).transpose())*D_ij.block(0, 2 * p, 2, 2)*AM_ij.block(0, p*u_M, 2, u_M);
				}
				NCC_j.block(0, 6 * j, 6, 6) = NCC_j.block(0, 6 * j, 6, 6) + (AC_ij.block(0, 6 * p, 2, 6).transpose())*D_ij.block(0, 2 * p, 2, 2)*AC_ij.block(0, 6 * p, 2, 6);
				UC_j.block(0, j, 6, 1) = UC_j.block(0, j, 6, 1) + (AC_ij.block(0, 6 * p, 2, 6).transpose())*D_ij.block(0, 2 * p, 2, 2)*W1_ij.block(0, p * 1, 2, 1);
			}
		}

		UC_j = -UC_j;

		MatrixXd NCC;
		NCC.setZero(u_o, u_o);

		for (int j = 0; j<n_cams; j++) {
			NCC.block(j * 6, j * 6, 6, 6) = NCC_j.block(0, 6 * j, 6, 6);
		}
		NCC_j.resize(0, 0); //this is to free the allocated memory space to a variable that we won't use anymore


		MatrixXd UC;
		UC.setZero(u_o, 1);
		for (int j = 0; j<n_cams; j++) {
			UC.block(6 * j, 0, 6, 1) = UC_j.block(0, j, 6, 1);
		}
		UC_j.resize(0, 0);

		MatrixXd NMG_i;
		if (u_M > 0) {
			NMG_i.setZero(u_M, 3 * n_points);//sigma_over_j(AM_ij'*D_ij*AG_ij)
		}
		MatrixXd NMGt_i;
		if (u_M > 0) {
			NMGt_i.setZero(3, u_M*n_points);//sigma_over_j(AG_ij'*D_ij*AM_ij)
		}

		MatrixXd NGGinv_i;
		NGGinv_i.setZero(3, 3 * n_points);//inv(sum_over_j(AG_ij'*P_ij*AG_ij))
		MatrixXd UG_i;
		UG_i.setZero(3, 1 * n_points);//sum_over_j(AG_ij'*P_ij*W1_ij)


		for (int i = 0; i<n_points; i++) {
			int strt = Adress_p_i(i);
			int nd = Adress_p_i(i + 1) - 1;

			for (int p = strt; p <= nd; p++) { //sum over j
				if (u_M > 0) {
					NMG_i.block(0, 3 * i, u_M, 3) = NMG_i.block(0, 3 * i, u_M, 3) + (AM_ij.block(0, p*u_M, 2, u_M).transpose())*D_ij.block(0, 2 * p, 2, 2)*AG_ij.block(0, 3 * p, 2, 3);
					NMGt_i.block(0, u_M*i, 3, u_M) = NMGt_i.block(0, u_M*i, 3, u_M) + (AG_ij.block(0, 3 * p, 2, 3).transpose())*D_ij.block(0, 2 * p, 2, 2)*AM_ij.block(0, p*u_M, 2, u_M);
				}
				NGGinv_i.block(0, 3 * i, 3, 3) = NGGinv_i.block(0, 3 * i, 3, 3) + (AG_ij.block(0, 3 * p, 2, 3).transpose())*D_ij.block(0, 2 * p, 2, 2)*AG_ij.block(0, 3 * p, 2, 3);
				UG_i.block(0, i, 3, 1) = UG_i.block(0, i, 3, 1) + (AG_ij.block(0, 3 * p, 2, 3).transpose())*D_ij.block(0, 2 * p, 2, 2)*W1_ij.block(0, p, 2, 1);
			}
		}

		for (int i = 0; i<n_points; i++) {
			NGGinv_i.block(0, 3 * i, 3, 3) = NGGinv_i.block(0, 3 * i, 3, 3).inverse();
		}

		UG_i = -UG_i;

		////////// For inner constraints ///////////////
		MatrixXd Mo_j;
		Mo_j.setZero(7, 6 * n_cams); //inner constraints with regard to EOPs of image j
		MatrixXd Mg_i;
		Mg_i.setZero(7, 3 * n_points); //inner constrains with regard to Gi

		for (int j = 0; j<n_cams; j++) {
			int sens_j = (int)camera_params.sensor_ID(j, 0);

			double Omega0 = Bundestim_ori(j * 5 + 4, 0);
			double Phi0 = Bundestim_ori(j * 5 + 4, 1);
			double Kappa0 = Bundestim_ori(j * 5 + 4, 2);

			double Xo0 = Bundestim_ori(j * 5 + 3, 0);
			double Yo0 = Bundestim_ori(j * 5 + 3, 1);
			double Zo0 = Bundestim_ori(j * 5 + 3, 2);


			Mo_j.block(0, 6 * j, 1, 6) << 0, 0, 0, 1, 0, 0;
			Mo_j.block(1, 6 * j, 1, 6) << 0, 0, 0, 0, 1, 0;
			Mo_j.block(2, 6 * j, 1, 6) << 0, 0, 0, 0, 0, 1;
			Mo_j.block(3, 6 * j, 1, 6) << 1, 0, 0, 0, -Zo0, Yo0;
			Mo_j.block(4, 6 * j, 1, 6) << sin(Omega0)*tan(Phi0), cos(Omega0), -sin(Omega0) / cos(Phi0), Zo0, 0, -Xo0;
			Mo_j.block(5, 6 * j, 1, 6) << -cos(Omega0)*tan(Phi0), sin(Omega0), cos(Omega0) / cos(Phi0), -Yo0, Xo0, 0;
			Mo_j.block(6, 6 * j, 1, 6) << 0, 0, 0, Xo0, Yo0, Zo0;//this is for scale; if you have a measured distance, then remove this and add the distance observation
		}

		for (int i = 0; i<n_points; i++) {

			Mg_i.block(0, 3 * i, 1, 3) << 1, 0, 0;
			Mg_i.block(1, 3 * i, 1, 3) << 0, 1, 0;
			Mg_i.block(2, 3 * i, 1, 3) << 0, 0, 1;
			Mg_i.block(3, 3 * i, 1, 3) << 0, -ground(i, 2), ground(i, 1);
			Mg_i.block(4, 3 * i, 1, 3) << ground(i, 2), 0, -ground(i, 0);
			Mg_i.block(5, 3 * i, 1, 3) << -ground(i, 1), ground(i, 0), 0;
			Mg_i.block(6, 3 * i, 1, 3) << ground(i, 0), ground(i, 1), ground(i, 2); //this is for scale; if you have a measured distance, then remove this and add the distance observation
		}


		////////////////////////////////////////////////
		MatrixXd KM;
		MatrixXd Zb_i;
		MatrixXd Zg_i;

		if (u_M > 0) {
			Zb_i.setZero(3, u_M*n_points); //HGinv_i*HCGt_i
			for (int i = 0; i < n_points; i++) {
				Zb_i.block(0, u_M*i, 3, u_M) = NGGinv_i.block(0, 3 * i, 3, 3)*NMGt_i.block(0, u_M*i, 3, u_M);
			}

			KM.setZero(u_M, u_M); //inv(KM)=inv(NMM-(Sum_over_i(NMG_i*Zb_i)))
			tmp.setZero(u_M, u_M);
			for (int i = 0; i < n_points; i++) {
				tmp = tmp + NMG_i.block(0, 3 * i, u_M, 3)*Zb_i.block(0, u_M*i, 3, u_M);
			}
			KM = (NMM - tmp).inverse();

			Zg_i.setZero(3, u_M*n_points);
			for (int i = 0; i < n_points; i++) {
				Zg_i.block(0, u_M*i, 3, u_M) = Zb_i.block(0, u_M*i, 3, u_M)*KM;
			}

			Zb_i.resize(0, 0);
		}

		MatrixXd KC_j;
		if (u_M > 0) {
			KC_j.setZero(u_M, 6 * n_cams); //Zh_j-NMC_j
			for (int j = 0; j < n_cams; j++) {
				int strt = Adress_i_p(j);
				int nd = Adress_i_p(j + 1) - 1;
				MatrixXd Zh_j = MatrixXd::Zero(u_M, 6);
				for (int indx = strt; indx <= nd; indx++) { //sigma ruye i
					int p = Index_i_p(indx);
					int i = Index_p_i(p, 0); //for i including GCPs, the right side is zero anyways because of Ag_ij(:,:,p)
					Zh_j = Zh_j + NMG_i.block(0, 3 * i, u_M, 3)*(NGGinv_i.block(0, 3 * i, 3, 3)*((AG_ij.block(0, 3 * p, 2, 3).transpose())*D_ij.block(0, 2 * p, 2, 2)*AC_ij.block(0, 6 * p, 2, 6)));
				}
				KC_j.block(0, 6 * j, u_M, 6) = Zh_j - NMC_j.block(0, j * 6, u_M, 6);
			}

			NMC_j.resize(0, 0);
		}
		////////////////

		MatrixXd KK;
		if (u_M > 0) {
			KK.setZero(u_M, 7);
			for (int i = 0; i < n_points; i++) {
				KK = KK + NMG_i.block(0, 3 * i, u_M, 3)*(NGGinv_i.block(0, 3 * i, 3, 3)*(Mg_i.block(0, 3 * i, 7, 3).transpose()));
			}
		}
		MatrixXd UK_i = MatrixXd::Zero(3, 7 * n_points);
		for (int i = 0; i < n_points; i++) {
			if (u_M > 0) {
				UK_i.block(0, 7 * i, 3, 7) = -Zg_i.block(0, u_M*i, 3, u_M)*KK - (NGGinv_i.block(0, 3 * i, 3, 3)*(Mg_i.block(0, 3 * i, 7, 3).transpose()));
			}
			else {
				UK_i.block(0, 7 * i, 3, 7) = -(NGGinv_i.block(0, 3 * i, 3, 3)*(Mg_i.block(0, 3 * i, 7, 3).transpose()));
			}
		}


		///////////////

		MatrixXd KN;
		if (u_M > 0) {
			KN.setZero(u_M, 1);
			tmp.setZero(u_M, 1);
			for (int i = 0; i < n_points; i++) {
				tmp = tmp + NMG_i.block(0, 3 * i, u_M, 3)*(NGGinv_i.block(0, 3 * i, 3, 3)*UG_i.block(0, i, 3, 1));
			}
			KN = UM - tmp;
		}

		MatrixXd HG_i = MatrixXd::Zero(3, 1 * n_points);
		for (int i = 0; i < n_points; i++) {
			if (u_M > 0) {
				HG_i.block(0, i, 3, 1) = (NGGinv_i.block(0, 3 * i, 3, 3)*UG_i.block(0, i, 3, 1)) - (Zg_i.block(0, u_M*i, 3, u_M)*KN);
			}
			else {
				HG_i.block(0, i, 3, 1) = (NGGinv_i.block(0, 3 * i, 3, 3)*UG_i.block(0, i, 3, 1));
			}
		}

		UG_i.resize(0, 0);
		if (u_M > 0) {
			NMGt_i.resize(0, 0);
			NMG_i.resize(0, 0);
			UM.resize(0, 0);
		}

		MatrixXd Zm = MatrixXd::Zero(7, 7);
		for (int i = 0; i < n_points; i++) {
			Zm = Zm + Mg_i.block(0, 3 * i, 7, 3)*UK_i.block(0, 7 * i, 3, 7);
		}
		Zm = -Zm.inverse();


		MatrixXd SN = MatrixXd::Zero(7, 1);
		tmp.setZero(7, 1);
		for (int i = 0; i < n_points; i++) {
			tmp = tmp + Mg_i.block(0, 3 * i, 7, 3)*HG_i.block(0, i, 3, 1);
		}
		SN = Zm * tmp;

		//cout << "reached above ZS" << endl;
		MatrixXd SO_j = MatrixXd::Zero(7, 6 * n_cams);
		MatrixXd ZS = MatrixXd::Zero(6 * n_cams, 6 * n_cams);
		MatrixXd HGC_ij;

		for (int j2 = 0; j2<n_cams; j2++) {


			HGC_ij.setZero(3, 6 * n_points);

			int strt2 = Adress_i_p(j2);
			int nd2 = Adress_i_p(j2 + 1) - 1;
			if (u_M > 0) {
				for (int i = 0; i < n_points; i++) {
					HGC_ij.block(0, 6 * i, 3, 6) = -Zg_i.block(0, u_M*i, 3, u_M)*KC_j.block(0, 6 * j2, u_M, 6);
				}
			}

			for (int indx = strt2; indx <= nd2; indx++) {
				int p2 = Index_i_p(indx);
				int i = Index_p_i(p2, 0);
				if (u_M > 0) {
					HGC_ij.block(0, 6 * i, 3, 6) = -Zg_i.block(0, u_M*i, 3, u_M)*KC_j.block(0, 6 * j2, u_M, 6) - (NGGinv_i.block(0, 3 * i, 3, 3)*((AG_ij.block(0, 3 * p2, 2, 3).transpose())*D_ij.block(0, 2 * p2, 2, 2)*AC_ij.block(0, 6 * p2, 2, 6)));

				}
				else {
					HGC_ij.block(0, 6 * i, 3, 6) = -(NGGinv_i.block(0, 3 * i, 3, 3)*((AG_ij.block(0, 3 * p2, 2, 3).transpose())*D_ij.block(0, 2 * p2, 2, 2)*AC_ij.block(0, 6 * p2, 2, 6)));

				}
			}
			tmp.setZero(7, 6);
			for (int i = 0; i < n_points; i++) {
				tmp = tmp + Mg_i.block(0, 3 * i, 7, 3)*HGC_ij.block(0, 6 * i, 3, 6);
			}
			SO_j.block(0, 6 * j2, 7, 6) = Zm * (tmp + Mo_j.block(0, 6 * j2, 7, 6));

			for (int j1 = 0; j1<n_cams; j1++) {
				int strt1 = Adress_i_p(j1);
				int nd1 = Adress_i_p(j1 + 1) - 1;
				tmp.setZero(6, 6);
				for (int indx = strt1; indx <= nd1; indx++) {
					int p1 = Index_i_p(indx);
					int i = Index_p_i(p1, 0);
					tmp = tmp + ((AC_ij.block(0, 6 * p1, 2, 6).transpose())*D_ij.block(0, 2 * p1, 2, 2)*AG_ij.block(0, 3 * p1, 2, 3))* (HGC_ij.block(0, 6 * i, 3, 6) + UK_i.block(0, 7 * i, 3, 7)*SO_j.block(0, 6 * j2, 7, 6));
				}
				ZS.block(j1 * 6, j2 * 6, 6, 6) = tmp;
			}
		}

		HGC_ij.resize(0, 0);
		Mg_i.resize(0, 0);
		Zm.resize(0, 0);

		MatrixXd MOt = Mo_j.transpose();
		Mo_j.resize(0, 0);

		MatrixXd ZR;
		MatrixXd ZO;
		if (u_M > 0) {
			MatrixXd Zq = KN + KK * SN;
			ZR.setZero(6 * n_cams, 1);
			for (int j = 0; j < n_cams; j++) {
				ZR.block(j * 6, 0, 6, 1) = (NMCt_j.block(0, u_M*j, 6, u_M)*KM)*Zq;
			}
			Zq.resize(0, 0);


			ZO.setZero(6 * n_cams, 6 * n_cams);
			for (int j1 = 0; j1<n_cams; j1++) {
				for (int j2 = 0; j2<n_cams; j2++) {
					ZO.block(j1 * 6, j2 * 6, 6, 6) = (NMCt_j.block(0, u_M*j1, 6, u_M)*KM)* (KC_j.block(0, 6 * j2, u_M, 6) + KK * SO_j.block(0, 6 * j2, 7, 6));
				}
			}
			NMCt_j.resize(0, 0);
		}

		MatrixXd ZU = MatrixXd::Zero(6 * n_cams, 1);
		for (int j = 0; j<n_cams; j++) {
			int strt = Adress_i_p(j);
			int nd = Adress_i_p(j + 1) - 1;
			tmp.setZero(6, 1);
			for (int indx = strt; indx <= nd; indx++) {
				int p = Index_i_p(indx);
				int i = Index_p_i(p, 0); //for i including GCPs the right side is zero 
				tmp = tmp + ((AC_ij.block(0, 6 * p, 2, 6).transpose())*D_ij.block(0, 2 * p, 2, 2)*AG_ij.block(0, 3 * p, 2, 3))*(HG_i.block(0, i, 3, 1) + UK_i.block(0, 7 * i, 3, 7)*SN);
			}
			ZU.block(j * 6, 0, 6, 1) = tmp;
		}

		MatrixXd SO = SO_j;
		SO_j.resize(0, 0);

		MatrixXd OC;
		MatrixXd OB;
		if (u_M > 0) {
			OC = NCC + (ZO + ZS + MOt * SO);
			OB = (UC - ZR - ZU - MOt * SN);
		}
		else {
			OC = NCC + (ZS + MOt * SO);
			OB = UC - ZU - MOt * SN;
		}

		ZU.resize(0, 0);
		NCC.resize(0, 0);
		UC.resize(0, 0);
		ZS.resize(0, 0);
		MOt.resize(0, 0);
		if (u_M > 0) {
			ZO.resize(0, 0);
			ZR.resize(0, 0);
		}
		FullPivHouseholderQR<MatrixXd> qr = OC.fullPivHouseholderQr();

		int qrrank = qr.rank();
		if (qrrank<(n_cams * 6)) {
			//throw invalid_argument( "Scale problem!" );
			Write_Mat("OC.txt", OC, 15);
			cout << "WARNING! rank of (OC)=" << qrrank << ", n_cams*6=" << n_cams * 6 << endl;
			system("pause");
			return;
		}


		MatrixXd deltaxcap_C = qr.solve(OB); //corrections to EOPs
		if (u_M == 0) {
			Sigma_Xcap_C.setZero(6 * n_cams, 6 * n_cams);
			Sigma_Xcap_C = qr.inverse();
		}

		OC.resize(0, 0);
		OB.resize(0, 0);

		MatrixXd deltatemp_K2 = SN + SO * deltaxcap_C; //Lagrange multipliers
		SN.resize(0, 0);
		SO.resize(0, 0);


		MatrixXd KC;
		if (u_M > 0) {
			KC.setZero(u_M, 6 * n_cams);
			for (int j = 0; j < n_cams; j++) {
				KC.block(0, j * 6, u_M, 6) = KC_j.block(0, 6 * j, u_M, 6);
			}
			KC_j.resize(0, 0);
		}

		MatrixXd deltaxcap_M;
		if (u_M > 0) {
			deltaxcap_M = KM * (KC*deltaxcap_C + KK * deltatemp_K2 + KN); //Corrections to IOPs

			KM.resize(0, 0);
			KN.resize(0, 0);
			KK.resize(0, 0);
		}

		MatrixXd deltaxcap_G = MatrixXd::Zero(n_points, 3); //Corrections to 3D coords
		MatrixXd tmp1;
		if (u_M == 0) {
			Sigma_Xcap_G.setZero(3, 3 * n_points);
		}

		for (int i = 0; i<n_points; i++) {

			//Calculate HGC_ij
			int strt = Adress_p_i(i);
			int nd = Adress_p_i(i + 1) - 1;
			MatrixXd HGC_ij = MatrixXd::Zero(3, 6 * n_cams);

			if (u_M > 0) {
				for (int j = 0; j < n_cams; j++) {
					HGC_ij.block(0, j * 6, 3, 6) = -Zg_i.block(0, u_M*i, 3, u_M)*KC.block(0, j * 6, u_M, 6);
				}
			}

			for (int p = strt; p <= nd; p++) {
				int j = Index_p_i(p, 1);
				if (u_M > 0) {
					HGC_ij.block(0, j * 6, 3, 6) = -Zg_i.block(0, u_M*i, 3, u_M)*KC.block(0, j * 6, u_M, 6) - (NGGinv_i.block(0, 3 * i, 3, 3)*((AG_ij.block(0, 3 * p, 2, 3).transpose())*D_ij.block(0, 2 * p, 2, 2)*AC_ij.block(0, 6 * p, 2, 6)));
				}
				else {
					HGC_ij.block(0, j * 6, 3, 6) = -(NGGinv_i.block(0, 3 * i, 3, 3)*((AG_ij.block(0, 3 * p, 2, 3).transpose())*D_ij.block(0, 2 * p, 2, 2)*AC_ij.block(0, 6 * p, 2, 6)));
				}
			}
			tmp1.setZero(3, 1);
			for (int j = 0; j<n_cams; j++) {
				tmp1 = tmp1 + HGC_ij.block(0, j * 6, 3, 6)*deltaxcap_C.block(j * 6, 0, 6, 1);
			}
			tmp = (HG_i.block(0, i, 3, 1) + tmp1 + UK_i.block(0, 7 * i, 3, 7)*deltatemp_K2);
			deltaxcap_G.row(i) = tmp.transpose();

			if (u_M == 0) {
				MatrixXd temp_s = MatrixXd::Zero(3, 3);
				for (int j = 0; j<n_cams; j++) {
					MatrixXd tmp2 = MatrixXd::Zero(3, 6);
					for (int j2 = 0; j2<n_cams; j2++) {
						tmp2 = tmp2 + HGC_ij.block(0, j2 * 6, 3, 6)*Sigma_Xcap_C.block(6 * j2, 6 * j, 6, 6);
					}
					temp_s = temp_s + tmp2 * (HGC_ij.block(0, j * 6, 3, 6).transpose())*NGGinv_i.block(0, 3 * i, 3, 3);
				}
				Sigma_Xcap_G.block(0, 3 * i, 3, 3) = NGGinv_i.block(0, 3 * i, 3, 3) + temp_s;
			}

			HGC_ij.resize(0, 0);
		}

		HG_i.resize(0, 0);
		NGGinv_i.resize(0, 0);
		if (u_M > 0) {
			KC.resize(0, 0);
		}
		UK_i.resize(0, 0);

		//update unknowns
		maxcorrection = 0.0;
		for (int i = 0; i<n_points; i++) {
			ground.row(i) = ground.row(i) + deltaxcap_G.row(i);
			double maxvaltemp = (deltaxcap_G.row(i).array().abs()).maxCoeff();
			if (maxvaltemp>maxcorrection) {
				maxcorrection = maxvaltemp;
			}
		}


		for (int j = 0; j<n_cams; j++) {
			Bundestim_ori.row(j * 5 + 4) = Bundestim_ori.row(j * 5 + 4) + deltaxcap_C.block(j * 6, 0, 3, 1).transpose(); //Omega Phi Kappa 
			Bundestim_ori.row(j * 5 + 3) = Bundestim_ori.row(j * 5 + 3) + deltaxcap_C.block(j * 6 + 3, 0, 3, 1).transpose(); //Xo Yo Zo
			Matrix3b3 r; //R_gtoi
			Rotation_g2i(Bundestim_ori(5 * j + 4, 0), Bundestim_ori(5 * j + 4, 1), Bundestim_ori(5 * j + 4, 2), r);
			Bundestim_ori.block(j * 5, 0, 3, 3) = r;
		}

		double maxvaltemp = (deltaxcap_C.array().abs()).maxCoeff();
		if (maxvaltemp>maxcorrection) {
			maxcorrection = maxvaltemp;
		}

		int sens_j = 0; 
		int places=0;
			
		K1_0(sens_j, 0) = K1_0(sens_j, 0) + deltaxcap_M(10 * places, 0);
		K2_0(sens_j, 0) = K2_0(sens_j, 0) + deltaxcap_M(10 * places + 1, 0);
		K3_0(sens_j, 0) = K3_0(sens_j, 0) + deltaxcap_M(10 * places + 2, 0);
		P1_0(sens_j, 0) = P1_0(sens_j, 0) + deltaxcap_M(10 * places + 3, 0);
		P2_0(sens_j, 0) = P2_0(sens_j, 0) + deltaxcap_M(10 * places + 4, 0);
		S1_0(sens_j, 0) = S1_0(sens_j, 0) + deltaxcap_M(10 * places + 5, 0);
		S2_0(sens_j, 0) = S2_0(sens_j, 0) + deltaxcap_M(10 * places + 6, 0);
		xpp_0(sens_j, 0) = xpp_0(sens_j, 0) + deltaxcap_M(10 * places + 7, 0);
		ypp_0(sens_j, 0) = ypp_0(sens_j, 0) + deltaxcap_M(10 * places + 8, 0);
		f_l_0(sens_j, 0) = f_l_0(sens_j, 0) + deltaxcap_M(10 * places + 9, 0);
	


		if (u_M > 0) {
			maxvaltemp = (deltaxcap_M.array().abs()).maxCoeff();
			if (maxvaltemp > maxcorrection) {
				maxcorrection = maxvaltemp;
			}
		}


		//Estimate the residual vectors and the variance factor
		Sigma_02_cap = 0.0;
		int maxindx = 0;
		double maxv = 0.0;
		for (int p = 0; p<n_eq; p++) {
			int i = Index_p_i(p, 0);
			int j = Index_p_i(p, 1);
			int sens_j = (int)camera_params.sensor_ID(j, 0);
			Matrix2b1 Zw_ij;
			if (u_M > 0) {
				Zw_ij = AM_ij.block(0, p*u_M, 2, u_M)*deltaxcap_M + AC_ij.block(0, 6 * p, 2, 6)*deltaxcap_C.block(j * 6, 0, 6, 1) + AG_ij.block(0, 3 * p, 2, 3)*(deltaxcap_G.row(i).transpose()) + W1_ij.block(0, p, 2, 1);
			}
			else {
				Zw_ij = AC_ij.block(0, 6 * p, 2, 6)*deltaxcap_C.block(j * 6, 0, 6, 1) + AG_ij.block(0, 3 * p, 2, 3)*(deltaxcap_G.row(i).transpose()) + W1_ij.block(0, p, 2, 1);
			}
			Matrix2b1 Vcap_ij = -Zy_ij.block(0, 2 * p, 2, 2)*(D_ij.block(0, 2 * p, 2, 2)*Zw_ij);


			if ((Vcap_ij.array().abs()).maxCoeff()>maxv) {
				maxindx = p;
				maxv = (Vcap_ij.array().abs()).maxCoeff();
			}

			MatrixXd weightmatrix;
			weightmatrix.setZero(2, 2);
			weightmatrix << 1 / (Sigmaxy_p(p, 0)*Sigmaxy_p(p, 0)), 0,
				0, 1 / (Sigmaxy_p(p, 0)*Sigmaxy_p(p, 0));//this is the p_ij
			Sigma_02_cap = Sigma_02_cap + Vcap_ij.transpose()*weightmatrix*Vcap_ij;
		}

		Sigma_02_cap = Sigma_02_cap / (n_eq * 2 + 7 - n_points * 3 - n_cams * 6 - u_M);

		double xdv = Obspixels_p(maxindx, 0);
		double ydv = Obspixels_p(maxindx, 1);
		cout << "		Largest errors happen at x,y: " << xdv << " , " << ydv << endl;
		

		cout << "     MaxV=" << maxv << endl;
		Sigma_Xcap_C = Sigma_02_cap * Sigma_Xcap_C;
		Sigma_Xcap_G = Sigma_02_cap * Sigma_Xcap_G;

		cout << "     Maximum parameter correction=" << maxcorrection << endl;
		cout << "     Sigma_02_cap=" << Sigma_02_cap << endl;

		deltaxcap_C.resize(0, 0);
		deltaxcap_G.resize(0, 0);
		W1_ij.resize(0, 0);
		AC_ij.resize(0, 0);
		AG_ij.resize(0, 0);
		D_ij.resize(0, 0);
		Zy_ij.resize(0, 0);
		if (u_M>0) {
			AM_ij.resize(0, 0);
			deltaxcap_M.resize(0, 0);
		}

	}//end of bundle loop

	for (int i = 0; i<n_cams; i++) {
		Bundestim_ori.row(i * 5 + 3) = Bundestim_ori.row(i * 5 + 3) + delta_trans;// Xo Yo Zo
	}

	//object-space coordinates of points
	ground.col(0) = ground.col(0) + delta_trans(0, 0)*(MatrixXd::Ones(n_points, 1));
	ground.col(1) = ground.col(1) + delta_trans(0, 1)*(MatrixXd::Ones(n_points, 1));
	ground.col(2) = ground.col(2) + delta_trans(0, 2)*(MatrixXd::Ones(n_points, 1));

	//returning back the IO parameters to the non-scaled position
	IO_scalfact = 1 / IO_scalfact;
	PS_0 = PS_0 * IO_scalfact;
	f_l_0 = f_l_0 * IO_scalfact;
	xpp_0 = xpp_0 * IO_scalfact;
	ypp_0 = ypp_0 * IO_scalfact;
	K1_0 = K1_0 / (pow(IO_scalfact, (int)2)); K2_0 = K2_0 / (pow(IO_scalfact, (int)4)); K3_0 = K3_0 / (pow(IO_scalfact, (int)6));
	P1_0 = P1_0 / IO_scalfact; P2_0 = P2_0 / IO_scalfact;

	//updating camera parameters with the new estimations of EOPs and IOPs
	camera_params.ResetEOPs(Bundestim_ori);
	int sens_j = 0;
	camera_params.ResetIOPs0(sens_j, true, PS_0(sens_j, 0), Cn_0(sens_j, 0), Rn_0(sens_j, 0), K1_0(sens_j, 0), K2_0(sens_j, 0), K3_0(sens_j, 0), P1_0(sens_j, 0), P2_0(sens_j, 0), S1_0(sens_j, 0), S2_0(sens_j, 0), xpp_0(sens_j, 0), ypp_0(sens_j, 0), f_l_0(sens_j, 0));
	

	cout << endl << "DoF=" << (n_eq * 2 + 7 - n_points * 3 - n_cams * 6 - u_M) << endl;

};
