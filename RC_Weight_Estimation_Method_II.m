function [W_wing, W_HT, W_VT, W_fuselage, W_structure, weight_details] = RC_Weight_Estimation_Method_II()
%% RC Aircraft Component Weight Estimation - Class II Method

%% Create a parameters structure
params = struct();
    
    % ----- MATERIAL PROPERTIES -----
params.rho_balsa = 5.5;        % Balsa wood (lb/ft^3) - typical range: 4-7
params.rho_basswood = 22.5;    % Basswood (lb/ft^3) - typical range: 20-25
params.rho_plywood = 35;       % Plywood (lb/ft^3)
params.rho_foam = 2.5;         % Foam for LE/TE (lb/ft^3)
params.rho_monokote = 0.8;     % Covering material (lb/ft^2)
params.rho_carbonfiber = 93.64;  %Carbon Fiber (lb/ft^3)
   
    % ----- WING SPECIFICATIONS -----
params.S_wing = 1.77;           % Wing area (ft^2)
params.b_wing = 3.0;           % Wing span (ft)
params.c_mean = 0.59;           % mean chord (ft)
params.tc_wing = 0.12;             % Thickness-to-chord ratio
    
    % Wing structure
params.n_ribs = 9;                    % Number of ribs
params.rib_thickness = 1/(16*12);     % Rib thickness (ft)
params.rib_material = params.rho_basswood;       % Rib material density
    
params.n_spars = 2;                    % Number of spars
params.spar_width = 0.25 / 12;         % Spar width (ft)
params.spar_height = 0.5 / 12;         % Spar height (ft)
params.spar_material = params.rho_carbonfiber;   % Spar material density
    
params.n_stringers = 8;                            % Number of stringers
params.stringer_area = (0.125 * 0.125) / 144;      % Cross-sectional area (ft^2)
params.stringer_material = params.rho_basswood;              % Stringer material density
    
params.skin_material = params.rho_monokote;   % Skin material density

% NO FOAM
params.use_LE_material = false;
params.LE_volume_ratio = 0;
params.LE_material = 0;
    
params.use_TE_material = false;
params.TE_volume_ratio = 0;
params.TE_material = 0;
    
    % ----- HORIZONTAL TAIL SPECIFICATIONS -----
params.S_Ht = 0.168;                     % HT area (ft^2)
params.b_HT = 14.4/12;                     % HT span (ft)
params.t_c_HT = 0.10;                  % Thickness-to-chord ratio
params.HT_structure_type = 'flat_plate';   % 'flat_plate' or 'wing_structure'
params.HT_plate_thickness = 0.288 / 12;    % Plate thickness if flat plate (ft)- NACA 0006
params.HT_material = params.rho_basswood;            % HT material density
    
    % ----- VERTICAL TAIL SPECIFICATIONS -----
params.S_VT = 0.168/2;                     % VT area (ft^2)
params.b_VT = 14.4/(12*2);                     % VT height (ft)
params.VT_plate_thickness = 0.288 / 12;    % Plate thickness (ft) - NACA 0006
params.VT_material = params.rho_basswood;            % VT material density
   
% ----- FUSELAGE SPECIFICATIONS -----
params.L_fuse = 17/12;           % Fuselage length (ft)
params.W_fuse = 3.125/12;           % Fuselage width (ft)
params.H_fuse = 3/12;           % Fuselage height (ft)
    
params.n_frames = 20/2;                   % Number of frames
params.frame_thickness = 1/ (16*12);   % Frame thickness (ft)
params.frame_material = params.rho_basswood;     % Frame material density
    
params.n_longeroons = 8;                        % Number of longerons/stringers
params.longeroon_area = (0.125 * 0.125) / 144;  % Cross-sectional area (ft^2)
params.longeroon_material = params.rho_basswood;       % Longeron material density
    
params.fuse_skin_material = params.rho_monokote;  % Fuselage skin material

%--------Booms----------
params.n_booms = 2;                             
params.boom_area = (3) / (8*144);
params.boom_material = params.rho_carbonfiber;


%     % ----- POWERPLANT -----
% params.W_motor = 0.5;          % Motor weight (lb)
% params.W_propeller = 0.2;      % Propeller weight (lb)
% params.W_ESC = 0.15;           % ESC weight (lb)
% params.W_battery = 1.2;        % Battery weight (lb)
% params.W_wires_power = 0.1;    % Power wires (lb)
% 
%     % ----- FLIGHT CONTROL SYSTEM -----
% params.W_receiver = 0.05;      % Receiver weight (lb)
% params.n_servos = 4;           % Number of servos
% params.W_servo = 0.08;         % Weight per servo (lb)
% params.W_wires_control = 0.05; % Control wires (lb)
% 
%     % ----- EXTERNAL EQUIPMENT -----
% params.W_landing_gear = 0.3;   % Landing gear (lb)
% params.W_extra_support = 0.2;  % Extra supports/brackets (lb)
% 
%     % ----- PAYLOAD -----
% params.W_payload = 2.0;        % Payload weight (lb)
    
%% Call calculation function with parameters structure
[W_wing, W_HT, W_VT, W_fuselage, W_structure, weight_details] = calculate_weights(params);
    
figure('Position', [100, 100, 900, 600]);
    
labels = {'Wing', 'H-Tail', 'V-Tail', 'Fuselage'};
weights = [W_wing, W_HT, W_VT, W_fuselage];
    
pie(weights, labels);
title('RC Aircraft Structural Weight Breakdown', 'FontSize', 14, 'FontWeight', 'bold');
colormap(jet(length(labels)));
end