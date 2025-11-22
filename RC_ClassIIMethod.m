% RC Class II Method
clc; close all; clear all;

%% Use Given Parameters to Estimate Weight
[W_wing, W_HT, W_VT, W_fuselage, W_booms, W_structure, W_avionics, W_powerplant, W_total, weight_details, cg_details] = RC_Weight_Estimation_Method_II();
fprintf('Total Structural Weight: %.4f lb\n', W_structure);
fprintf('Wing Weight: %.4f lb\n', W_wing);
fprintf('HT Weight: %.4f lb\n', W_HT);
fprintf('VT Weight: %.4f lb\n', W_VT);
fprintf('Fuselage Weight: %.4f lb\n', W_fuselage);
fprintf('Boom Weight: %.4f lb\n', W_booms);



