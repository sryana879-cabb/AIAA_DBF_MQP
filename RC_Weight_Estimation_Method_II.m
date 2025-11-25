function [W_wing, W_HT, W_VT, W_fuselage, W_booms, W_structure, W_avionics, W_powerplant, W_total, weight_details, cg_details] = RC_Weight_Estimation_Method_II()
%% RC Aircraft Component Weight Estimation - Class II Method with CG

%% Create a parameters structure
params = struct();
    
    % ----- MATERIAL PROPERTIES -----
params.rho_balsa = 5.5;        % Balsa wood (lb/ft^3) - typical range: 4-7
params.rho_basswood = 22.5;    % Basswood (lb/ft^3) - typical range: 20-25
params.rho_plywood = 35;       % Plywood (lb/ft^3)
params.rho_foam = 2.5;         % Foam for LE/TE (lb/ft^3)
params.rho_monokote = 0.0125;     % Covering material (lb/ft^2)
params.rho_carbonfiber = 93.64;  %Carbon Fiber (lb/ft^3)
   
    % ----- WING SPECIFICATIONS -----
params.S_wing = 2.557014;           % Wing area (ft^2) og = 1.77
params.b_wing = 3.198133;           % Wing span (ft)
params.c_mean = 0.79953;           % mean chord (ft), og = 0.59
params.tc_wing = 0.12;             % Thickness-to-chord ratio
    
    % Wing structure
params.n_ribs = 9;                    % Number of ribs
params.rib_thickness = 1/(16*12);     % Rib thickness (ft)
params.rib_material = params.rho_balsa;       % Rib material density
    
params.n_spars = 2;                    % Number of spars
params.spar_width = 0.25 / 12;         % Spar width (ft)
params.spar_height = 0.5 / 12;         % Spar height (ft)
params.spar_material = params.rho_basswood;   % Spar material density
    
params.n_stringers = 8;                            % Number of stringers
params.stringer_area = (0.125 * 0.125) / 144;      % Cross-sectional area (ft^2)
params.stringer_material = params.rho_basswood;     % Stringer material density
    
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
params.HT_plate_thickness = 0.27 / 12;    % Plate thickness if flat plate (ft)- NACA 0006
params.HT_material = params.rho_basswood;            % HT material density

    
    % ----- VERTICAL TAIL SPECIFICATIONS (for both tails)-----
params.S_VT = 49.5/144;                     % VT area (ft^2)
params.b_VT = 11/(12);                     % VT height (ft)
params.VT_plate_thickness = 0.27 / 12;    % Plate thickness (ft) - NACA 0006
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
params.L_boom = params.L_fuse;
params.boom_area = (3/8)^2/144;
params.boom_material = params.rho_carbonfiber;

%% Component Positions (From Wing LE)
% Vertical Tail
params.x_VT_LE = params.c_mean + params.L_boom;          % ft aft of wing LE
params.z_VT = params.b_VT*0.5;            % ft above wing LE height

% Horizontal Tail position
params.x_HT_LE = params.c_mean + params.L_boom;          % ft aft of wing LE
params.z_HT = 0;            % ft above wing LE height

% Fuselage position
params.x_fuse_nose = -0.25*params.L_fuse;     % ft, assumed 1/4 to get to nose
params.z_fuse_bottom = -params.H_fuse;   % ft (negative = below wing LE)

%Boom Position
params.x_boom_start = params.c_mean;    % ft (where booms start, likely at nose)
params.z_boom = 0.0;           % ft (boom height, likely near wing height)

   % ----- POWERPLANT POSITIONS -----
params.W_motor = 0.487222;
params.x_motor = 0;  % Forward of wing LE
params.z_motor = 0.0;

params.W_propeller = 0.075001;
params.x_propeller = 0;  % Even more forward
params.z_propeller = 0.0;

params.W_ESC = 0.401241;
params.x_ESC = 0.05;  % Just ahead of wing
params.z_ESC = 0;

params.W_battery = 0.529;
params.x_battery = 0.25 * params.c_mean;  % assumption
params.z_battery = 0;

    % ----- FLIGHT CONTROL SYSTEM POSITIONS -----
params.W_receiver = 0.035274;
params.x_receiver = 0.25 * params.c_mean;  % 25% chord
params.z_receiver = 0.0;

params.n_servos = 6;
params.W_servo = 0.0412;

%2 servo locations
c_HT = params.S_Ht/params.b_HT;
x_wing_servos = 2*(params.c_mean - 2.2/12) + 2*(params.c_mean-1.5/12);
x_tail_servos = params.x_HT_LE + 0.5*c_HT;

params.x_servos = (2*x_tail_servos + 4*x_wing_servos)/params.n_servos ; 
params.z_servos = 0;

    % ----- EXTERNAL EQUIPMENT -----
params.W_landing_gear = 0.3;
params.x_landing_gear = 0.20 * params.c_mean;
params.z_landing_gear = -params.H_fuse - 0.05;

params.W_extra_support = 0.2;

% ----- PAYLOAD -----
params.W_payload = 0.04375*3+0.375;
params.x_payload = 0.30 * params.c_mean;  
params.z_payload = 0.0;
    
%% Call calculate_weights
[W_wing, W_HT, W_VT, W_fuselage, W_booms, W_structure, W_avionics, W_powerplant, W_total, weight_details, cg_details] = calculate_weights(params);

% Check if cg_details was returned (for compatibility)
if ~exist('cg_details', 'var')
    error('CG details not calculated. Make sure calculate_weights_copy returns cg_details.');
end
    
%% Visualizations
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Weight breakdown pie chart
subplot(1, 3, 1);
labels = {'Wing', 'H-Tail', 'V-Tail', 'Fuselage', 'Booms', 'Powerplant', 'Avionics', 'Payload'};
weights = [W_wing, W_HT, W_VT, W_fuselage, W_booms, W_powerplant, W_avionics, params.W_payload];
pie(weights, labels);
title('RC Aircraft Structural Weight Breakdown', 'FontSize', 14, 'FontWeight', 'bold');
colormap(jet(length(labels)));

% Subplot 2: Structural Weight
subplot(1, 3, 2);
axis off;

% Create table data 
component_col = {'Wing'; 'H-Tail'; 'V-Tail'; 'Fuselage'; 'Booms'; 
                 'Powerplant'; 'Avionics'; ''; 'EMPTY'};
weight_col = {sprintf('%.3f', W_wing); sprintf('%.3f', W_HT); 
              sprintf('%.3f', W_VT); sprintf('%.3f', W_fuselage); 
              sprintf('%.3f', W_booms); sprintf('%.3f', W_powerplant);
              sprintf('%.3f', W_avionics); 
              ''; sprintf('%.3f', cg_details.aircraft_empty.weight)};
              
% Calculate weighted average CG for powerplant
x_cg_powerplant = (2*params.W_motor*params.x_motor + ...
                   2*params.W_propeller*params.x_propeller + ...
                   2*params.W_ESC*params.x_ESC + ...
                   params.W_battery*params.x_battery) / W_powerplant;

% Calculate weighted average CG for avionics
x_cg_avionics = (params.W_receiver*params.x_receiver + ...
                 params.n_servos*params.W_servo*params.x_servos) / W_avionics;

x_cg_col = {sprintf('%.1f', cg_details.wing.x_cg*12); 
            sprintf('%.1f', cg_details.HT.x_cg*12);
            sprintf('%.1f', cg_details.VT.x_cg*12); 
            sprintf('%.1f', cg_details.fuselage.x_cg*12);
            sprintf('%.1f', cg_details.booms.x_cg*12);
            sprintf('%.1f', x_cg_powerplant*12);
            sprintf('%.1f', x_cg_avionics*12);
            '';
            sprintf('%.1f', cg_details.aircraft_empty.x_cg*12)};

% Display as table
table_data = [component_col, weight_col, x_cg_col];
col_names = {'Component', 'Weight (lb)', 'x-CG (in)'};

% Create text table
y_start = 0.95;
y_step = 0.09; 

% Headers
text(0.1, y_start, col_names{1}, 'FontWeight', 'bold', 'FontSize', 10);
text(0.5, y_start, col_names{2}, 'FontWeight', 'bold', 'FontSize', 10);
text(0.75, y_start, col_names{3}, 'FontWeight', 'bold', 'FontSize', 10);
line([0.05 0.95], [y_start-0.02 y_start-0.02], 'Color', 'k', 'LineWidth', 2);

for i = 1:9
    y_pos = y_start - (i+0.3)*y_step;
    if i == 9
        line([0.05 0.95], [y_pos+0.03 y_pos+0.03], 'Color', 'k', 'LineWidth', 1);
        text(0.1, y_pos, table_data{i,1}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.5, y_pos, table_data{i,2}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.75, y_pos, table_data{i,3}, 'FontSize', 9, 'FontWeight', 'bold');
    elseif i == 8
        continue;
    else
        text(0.1, y_pos, table_data{i,1}, 'FontSize', 9);
        text(0.5, y_pos, table_data{i,2}, 'FontSize', 9);
        text(0.75, y_pos, table_data{i,3}, 'FontSize', 9);
    end
end
title('Empty Weight CG', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 3: Takeoff weight comparison
subplot(1, 3, 3);
axis off;

fprintf('\nCG SHIFT WITH PAYLOAD:\n');
fprintf('Empty:   %.2f in (%.1f%% MAC)\n', ...
    cg_details.aircraft_empty.x_cg*12, ...
    (cg_details.aircraft_empty.x_cg / params.c_mean) * 100);
fprintf('Takeoff: %.2f in (%.1f%% MAC)\n', ...
    cg_details.aircraft_takeoff.x_cg*12, ...
    (cg_details.aircraft_takeoff.x_cg / params.c_mean) * 100);
fprintf('Shift:   %.2f in (%.1f%% MAC)\n', ...
    (cg_details.aircraft_takeoff.x_cg - cg_details.aircraft_empty.x_cg)*12, ...
    ((cg_details.aircraft_takeoff.x_cg - cg_details.aircraft_empty.x_cg) / params.c_mean) * 100);

text(0.1, 0.8, 'CG COMPARISON', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.65, sprintf('Empty Weight: %.2f lb', cg_details.aircraft_empty.weight), 'FontSize', 10);
text(0.1, 0.55, sprintf('x-CG: %.2f in (%.1f%% MAC)', ...
    cg_details.aircraft_empty.x_cg*12, (cg_details.aircraft_empty.x_cg / params.c_mean) * 100), 'FontSize', 10);

text(0.1, 0.35, sprintf('Takeoff Weight: %.2f lb', cg_details.aircraft_takeoff.weight), 'FontSize', 10);
text(0.1, 0.25, sprintf('x-CG: %.2f in (%.1f%% MAC)', ...
    cg_details.aircraft_takeoff.x_cg*12, (cg_details.aircraft_takeoff.x_cg / params.c_mean) * 100), 'FontSize', 10);

text(0.1, 0.05, sprintf('CG Shift: %.2f in (%.1f%% MAC)', ...
    (cg_details.aircraft_takeoff.x_cg - cg_details.aircraft_empty.x_cg)*12, ...
    ((cg_details.aircraft_takeoff.x_cg - cg_details.aircraft_empty.x_cg) / params.c_mean) * 100), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');

end