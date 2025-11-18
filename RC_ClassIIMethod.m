% RC Class II Method
clc; close all; clear all;

%% Use Given Parameters to Estimate Weight
[W_wing, W_HT, W_VT, W_fuselage, W_structure, weight_details] = RC_Weight_Estimation_Method_II();
fprintf('Total Structural Weight: %.4f lb\n', W_structure);
fprintf('Total Wing Weight: %.4f lb\n', W_wing);
fprintf('Total HT Weight: %.4f lb\n', W_HT);
fprintf('Total VT Weight: %.4f lb\n', W_VT);
fprintf('Total Fuselage Weight: %.4f lb\n', W_fuselage);



