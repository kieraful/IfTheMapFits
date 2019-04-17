%ENGO500 If The Map Fits
%Bundle Adjustment - Derivatives for design matrices A and B
clear 
clc
%x
syms x_Sjb y_Sjb z_Sjb w_Sb phi_Sb K_Sb             double real  %boresight parameters, lidar to gnss/ins
syms n_xpg n_ypg n_zpg d_p                          double real   %components of a wall's normal vector, distance from lidar to point on wall
syms x_bjg y_bjg z_bjg w_bjg phi_bjg K_bjg          double real   %scan in ground ref frame

%l
syms x_sj y_sj z_sj                                 double real   %point coordinates OBSERVED with lidar
syms x_GPSj y_GPSj z_GPSj w_INSj phi_INSj K_INSj    double real   %OBSERVED GNSS ground coordinates and INS rotations
syms x_bs y_bs z_bs w_bs phi_bs K_bs                double real   % OBSERVED BORESIGHT OFFSETS, obs of x_Sjb y_Sjb z_Sjb


%other
syms A1
syms A2
syms A3
syms B1
syms B2
syms B3

%% Observation Equation Type 1

ObsEqn1 = [x_bjg; y_bjg; z_bjg; w_bjg; phi_bjg; K_bjg] ...
    - [x_GPSj; y_GPSj; z_GPSj; w_INSj; phi_INSj; K_INSj];
%Actually 6 equations, per scene.

A1 = simplify(jacobian(ObsEqn1, ...
    [x_Sjb y_Sjb z_Sjb w_Sb phi_Sb K_Sb ...
    n_xpg n_ypg n_zpg ...
    d_p ...
    x_bjg y_bjg z_bjg w_bjg phi_bjg K_bjg]));
% A1diff = diff(ObsEqn1(1,1), ...
%     x_Sjb, y_Sjb, z_Sjb, w_Sb, phi_Sb, K_Sb, ...
%     n_xpg, n_ypg, n_zpg, ...
%     d_p, ...
%     x_bjg, y_bjg, z_bjg, w_bjg, phi_bjg, K_bjg);
%6 by 16
A1args = argnames(A1);

B1 = simplify(jacobian(ObsEqn1, ...
    [x_sj y_sj z_sj ...
    x_GPSj y_GPSj z_GPSj w_INSj phi_INSj K_INSj]));
%6 by 9


%% Observation Equation Type 2

ObsEqn2 = n_xpg^2 + n_ypg^2 + n_zpg^2 - 1;

A2 = simplify(jacobian(ObsEqn2, ...
    [x_Sjb y_Sjb z_Sjb w_Sb phi_Sb K_Sb ...
    n_xpg n_ypg n_zpg ...
    d_p ...
    x_bjg y_bjg z_bjg w_bjg phi_bjg K_bjg]));
%1 by 16
A2args = argnames(A2);

B2 = simplify(jacobian(ObsEqn2, ...
    [x_sj y_sj z_sj ...
    x_GPSj y_GPSj z_GPSj w_INSj phi_INSj K_INSj]));
%1 by 9


%% Observation Equation Type 3
% Rx_bjg = [1 0 0; 0 cos(w_bjg) sin(w_bjg); 0 -sin(w_bjg) cos(w_bjg)]';
% Ry_bjg = [cos(phi_bjg) 0 -sin(phi_bjg); 0 1 0; sin(phi_bjg) 0 cos(phi_bjg)]';
% Rz_bjg = [cos(K_bjg) sin(K_bjg) 0; -sin(K_bjg) cos(K_bjg) 0; 0 0 1]';
% Rbjg=Rz_bjg*Ry_bjg*Rx_bjg;

%formula is for object to image. put in object to image angles (ground to
%IMU and IMU to  scan) and then transpose to get image to object

%with transpose
%without transpose
%with transpose, negated angles
%without transpose, negated angles
%no transpose, no negating

% Rx_bjg = [1 0 0; 0 cos(w_bjg) sin(w_bjg); 0 -sin(w_bjg) cos(w_bjg)];
% Ry_bjg = [cos(phi_bjg) 0 -sin(phi_bjg); 0 1 0; sin(phi_bjg) 0 cos(phi_bjg)];
% Rz_bjg = [cos(K_bjg) sin(K_bjg) 0; -sin(K_bjg) cos(K_bjg) 0; 0 0 1];
% Rbjg=Rz_bjg*Ry_bjg*Rx_bjg;
% %Rbjg=Rbjg';
% 
% Rx_Sb = [1 0 0; 0 cos(w_Sb) sin(w_Sb); 0 -sin(w_Sb) cos(w_Sb)];
% Ry_Sb = [cos(phi_Sb) 0 -sin(phi_Sb); 0 1 0; sin(phi_Sb) 0 cos(phi_Sb)];
% Rz_Sb = [cos(K_Sb) sin(K_Sb) 0; -sin(K_Sb) cos(K_Sb) 0; 0 0 1];
% RSb=Rz_Sb*Ry_Sb*Rx_Sb;
%RSb=RSb';

% Mimg2obj =
%  
% [                  cosK*cosphi,                 -cosphi*sinK,       sinphi]
% [ cosw*sinK + cosK*sinw*sinphi, cosK*cosw - sinK*sinw*sinphi, -cosphi*sinw]
% [ sinK*sinw - cosK*cosw*sinphi, cosK*sinw + cosw*sinK*sinphi,  cosw*cosphi]

% Rbjg = [cos(K_bjg)*cos(phi_bjg), -cos(phi_bjg)*sin(K_bjg),sin(phi_bjg);
%     cos(w_bjg)*sin(K_bjg) + cos(K_bjg)*sin(w_bjg)*sin(phi_bjg), cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(w_bjg)*sin(phi_bjg), -cos(phi_bjg)*sin(w_bjg);
%     sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg), cos(K_bjg)*sin(w_bjg) + cos(w_bjg)*sin(K_bjg)*sin(phi_bjg),  cos(w_bjg)*cos(phi_bjg)]';
% 
% RSb = [cos(K_Sb)*cos(phi_Sb), -cos(phi_Sb)*sin(K_Sb),sin(phi_Sb);
%     cos(w_Sb)*sin(K_Sb) + cos(K_Sb)*sin(w_Sb)*sin(phi_Sb), cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(w_Sb)*sin(phi_Sb), -cos(phi_Sb)*sin(w_Sb);
%     sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb), cos(K_Sb)*sin(w_Sb) + cos(w_Sb)*sin(K_Sb)*sin(phi_Sb),  cos(w_Sb)*cos(phi_Sb)]';

% Rbjg = [cos(-K_bjg)*cos(-phi_bjg), -cos(-phi_bjg)*sin(-K_bjg),sin(-phi_bjg);
%     cos(-w_bjg)*sin(-K_bjg) + cos(-K_bjg)*sin(-w_bjg)*sin(-phi_bjg), cos(-K_bjg)*cos(-w_bjg) - sin(-K_bjg)*sin(-w_bjg)*sin(-phi_bjg), -cos(-phi_bjg)*sin(-w_bjg);
%     sin(-K_bjg)*sin(-w_bjg) - cos(-K_bjg)*cos(-w_bjg)*sin(-phi_bjg), cos(-K_bjg)*sin(-w_bjg) + cos(-w_bjg)*sin(-K_bjg)*sin(-phi_bjg),  cos(-w_bjg)*cos(-phi_bjg)]';
% 
% RSb = [cos(-K_Sb)*cos(-phi_Sb), -cos(-phi_Sb)*sin(-K_Sb),sin(-phi_Sb);
%     cos(-w_Sb)*sin(-K_Sb) + cos(-K_Sb)*sin(-w_Sb)*sin(-phi_Sb), cos(-K_Sb)*cos(-w_Sb) - sin(-K_Sb)*sin(-w_Sb)*sin(-phi_Sb), -cos(-phi_Sb)*sin(-w_Sb);
%     sin(-K_Sb)*sin(-w_Sb) - cos(-K_Sb)*cos(-w_Sb)*sin(-phi_Sb), cos(-K_Sb)*sin(-w_Sb) + cos(-w_Sb)*sin(-K_Sb)*sin(-phi_Sb),  cos(-w_Sb)*cos(-phi_Sb)]';

Rbjg = [cos(K_bjg)*cos(phi_bjg), -cos(phi_bjg)*sin(K_bjg),sin(phi_bjg);
    cos(w_bjg)*sin(K_bjg) + cos(K_bjg)*sin(w_bjg)*sin(phi_bjg), cos(K_bjg)*cos(w_bjg) - sin(K_bjg)*sin(w_bjg)*sin(phi_bjg), -cos(phi_bjg)*sin(w_bjg);
    sin(K_bjg)*sin(w_bjg) - cos(K_bjg)*cos(w_bjg)*sin(phi_bjg), cos(K_bjg)*sin(w_bjg) + cos(w_bjg)*sin(K_bjg)*sin(phi_bjg),  cos(w_bjg)*cos(phi_bjg)];

RSb = [cos(K_Sb)*cos(phi_Sb), -cos(phi_Sb)*sin(K_Sb),sin(phi_Sb);
    cos(w_Sb)*sin(K_Sb) + cos(K_Sb)*sin(w_Sb)*sin(phi_Sb), cos(K_Sb)*cos(w_Sb) - sin(K_Sb)*sin(w_Sb)*sin(phi_Sb), -cos(phi_Sb)*sin(w_Sb);
    sin(K_Sb)*sin(w_Sb) - cos(K_Sb)*cos(w_Sb)*sin(phi_Sb), cos(K_Sb)*sin(w_Sb) + cos(w_Sb)*sin(K_Sb)*sin(phi_Sb),  cos(w_Sb)*cos(phi_Sb)];




ObsEqn3 = ( [x_bjg; y_bjg; z_bjg] ...
    + Rbjg * [x_Sjb; y_Sjb; z_Sjb] ...
    + Rbjg * RSb * [x_sj; y_sj; z_sj])' ...
    * [n_xpg; n_ypg; n_zpg] + d_p;

A3 = (jacobian(ObsEqn3, ...
    [x_Sjb y_Sjb z_Sjb w_Sb phi_Sb K_Sb ...
    n_xpg n_ypg n_zpg ...
    d_p ...
    x_bjg y_bjg z_bjg w_bjg phi_bjg K_bjg]));

%A3 is 1 by 16. For each equation that is type 1, there will be a row
%similar to this 1 by 16 in A (but with additional columns of zeros).
A3args = argnames(A3);

B3 = simplify(jacobian(ObsEqn3, ...
    [x_sj y_sj z_sj ...
    x_GPSj y_GPSj z_GPSj w_INSj phi_INSj K_INSj]));
%1 by 9


%% Simplifying ObsEqn3
syms Rbjg11 Rbjg12 Rbjg13 Rbjg21 Rbjg22 Rbjg23 Rbjg31 Rbjg32 Rbjg33 double real
syms RSb11 RSb12 RSb13 RSb21 RSb22 RSb23 RSb31 RSb32 RSb33 double real

Rbjg_sym = [Rbjg11, Rbjg12, Rbjg13; Rbjg21, Rbjg22, Rbjg23; Rbjg31, Rbjg32, Rbjg33];
RSb_sym = [RSb11, RSb12, RSb13; RSb21, RSb22, RSb23; RSb31, RSb32, RSb33];

ObsEqn3sym = simplify(( [x_bjg; y_bjg; z_bjg] ...
    + Rbjg_sym * [x_Sjb; y_Sjb; z_Sjb] ...
    + Rbjg_sym * RSb_sym * [x_sj; y_sj; z_sj])' ...
    * [n_xpg; n_ypg; n_zpg] + d_p);

% A3sym = simplify(jacobian(ObsEqn3sym, ...
%     [x_Sjb y_Sjb z_Sjb w_Sb phi_Sb K_Sb ...
%     n_xpg n_ypg n_zpg ...
%     d_p ...
%     x_bjg y_bjg z_bjg w_bjg phi_bjg K_bjg]));


% For bjg
% Deriv of ObsEqn3 wrt elements of Rbjg rotation matrix
D_Eqn3_Rbjg = simplify(jacobian(ObsEqn3sym,[Rbjg11, Rbjg12, Rbjg13, Rbjg21, Rbjg22, Rbjg23, Rbjg31, Rbjg32, Rbjg33]));
%Deriv of Rbjg elements wrt w_bjg, phi_bjg, K_bjg
D_Rbjg = simplify(jacobian(reshape(Rbjg,[9,1]), [w_bjg,phi_bjg,K_bjg]));
A3sym = D_Eqn3_Rbjg*D_Rbjg;


% For Sb
% Deriv of ObsEqn3 wrt elements of RSb rotation matrix
D_Eqn3_RSb = simplify(jacobian(ObsEqn3sym,[RSb11, RSb12, RSb13, RSb21, RSb22, RSb23, RSb31, RSb32, RSb33]));
%Deriv of RSb elements wrt w_Sb, phi_Sb, K_Sb
D_RSb = simplify(jacobian(reshape(RSb,[9,1]), [w_Sb,phi_Sb,K_Sb]));
A3sym = D_Eqn3_RSb*D_RSb;


A3check=A3sym-A3(1,4:6);


%% Boresight offsets

ObsEqn4 = [x_Sjb; y_Sjb; z_Sjb; w_Sb; phi_Sb; K_Sb] ...
    - [x_bs; y_bs; z_bs; w_bs; phi_bs; K_bs ];

A4 = (jacobian(ObsEqn4, ...
    [x_Sjb y_Sjb z_Sjb w_Sb phi_Sb K_Sb ...
    n_xpg n_ypg n_zpg ...
    d_p ...
    x_bjg y_bjg z_bjg w_bjg phi_bjg K_bjg]));
