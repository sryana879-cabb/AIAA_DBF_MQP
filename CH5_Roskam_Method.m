% Roskam method for CH5 
% Zero lift drag calculation

close all; clear all; clc; 

cw_mean = 10.4; % inches, mean chord of wing
cht_mean = 4; % inches, mean chord of HT
cvt_mean = 4.5; % inches, mean chord of VT
L_f = 17; % fuselage length
Re = 0.369e6; % Reynolds number for wing
Re_F = Re * L_f / cw_mean; % reynolds number of fuselage
Re_HT = Re * cht_mean / cw_mean; % reynolds number of HT
Re_VT = Re * cvt_mean / cw_mean; % reynolds number of VT
k = 0.00083; % surface roughness of cast iron

S = 3; % feet^2, wing planform area
S_wetW = S*2; % feet, wetted area of wing
S_wetHT = 80*2 / 144; % feet, wetted area of horizontal tail
S_wetVT = 20*2 / 144; % feet, wetted area of vertical tail
S_wetf = 198.9851 / 144; % feet, wetted area of fuselage
S_wetg = 0.22; % average wetted area of landing gear

% thickness to chord ratio t/c
tc_W = 1.248/10.4; 
tc_HT = 0.06; 
tc_VT = 0.06; 

% wing/fuselage interferenace factor from Fig 4.1
R_WF = 1.1; 
R_HTF = 1;
R_VTF = 1;
R_FF = 1.1;

% lift surface correction factor from Fig. 4.2
R_LSW = 1;
R_LSHT = 1;

% turbulent flate plate friction coefficient of the wing from Fig 4.3, increased 20% due to surface roughness
C_fW = 0.0055 * 1.2; 
C_fHT = 0.006 * 1.2;
C_fVT = 0.006 * 1.2;
C_fF = 0.0042 * 1.2;
Cf_g = 1.0; 
L = 2; % airfoil thickness location parameter from Fig. 4.4

% zero lift drag coefficient
C_D0W = R_WF*C_fW*S_wetW/S; % wing
C_D0HT = R_HTF*C_fHT*S_wetHT/S; % HT
C_D0VT = R_VTF*C_fVT*S_wetVT/S; % VT
C_D0F = R_FF*C_fF*S_wetf/S; % fuselage
C_D0G = 0.00446; % landing gear, estimated

% add coefficients of components
CD0_components = C_D0W + C_D0HT + C_D0VT + C_D0F + C_D0G;

% add 10% for CD_misc
CD_misc = 0.10 * CD0_components;       % 10% of all parasite drag
CD0_total = CD0_components + CD_misc;

fprintf("----- ZERO-LIFT DRAG RESULTS -----\n");
fprintf("CD0_wing       = %.5f\n", C_D0W);
fprintf("CD0_HT (flat)  = %.5f\n", C_D0HT);
fprintf("CD0_VT (flat)  = %.5f\n", C_D0VT);
fprintf("CD0_fuselage   = %.5f\n", C_D0F);
fprintf("CD0_landinggear = %.5f\n", C_D0G);
fprintf("TOTAL CD0      = %.5f\n", CD0_total);

