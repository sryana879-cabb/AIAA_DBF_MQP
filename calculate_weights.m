 function [W_wing, W_HT, W_VT, W_fuselage, W_structure, weight_details] = calculate_weights(p)
%% Wing Calculations
AR_wing = p.b_wing^2/p.S_wing;

%Ribs (with lightening factor for cutouts)
rib_chord_avg = p.c_mean;
rib_height = p.tc_wing*rib_chord_avg;
rib_perimeter = 2 *(rib_chord_avg+rib_height);
rib_area = rib_perimeter * p.rib_thickness;
W_ribs = p.n_ribs * rib_area * p.rib_material * 0.4;

%Spars
spar_length = p.b_wing;
spar_volume = p.n_spars * spar_length * p.spar_width * p.spar_height;
W_spars = spar_volume * p.spar_material;

%Stringers
stringer_length = p.b_wing;
W_stringers = p.n_stringers * p.stringer_area * stringer_length * p.stringer_material;

%Skin
W_skin = p.S_wing * p.skin_material;
wing_volume = p.S_wing * p.c_mean * p.tc_wing;

W_wing = W_ribs + W_spars + W_stringers + W_skin;

%% Horizontal Tail Calculations
c_HT = p.S_Ht/p.b_HT;

if strcmp(p.HT_structure_type, 'flat_plate')
    W_HT = p.S_Ht * p.HT_plate_thickness * p.HT_material;
else
    W_HT = p.S_Ht * 0.15;  % Simplified scaling for wing structure
end

%% Vertical Tail Calculations
c_VT = p.S_VT / p.b_VT;
W_VT = p.S_VT * p.VT_plate_thickness * p.VT_material;

%% Fuselage Calculations
% Frames (with lightening factor for cutouts)
frame_perimeter = 2*(p.W_fuse + p.H_fuse);
W_frames = p.n_frames * frame_perimeter * p.frame_thickness * p.frame_material * 0.4;

W_longeroons = p.n_longeroons * p.L_fuse * p.longeroon_area * p.longeroon_material;

fuse_surface_area = 2*(p.W_fuse * p.L_fuse + p.H_fuse * p.L_fuse);

W_fuse_skin = fuse_surface_area * p.fuse_skin_material;

W_booms = p.n_booms * p.L_fuse * p.boom_area * p.boom_material;
W_fuselage = W_frames + W_longeroons + W_fuse_skin + W_booms;

% %% Powerplant
% 
% W_powerplant = p.W_motor + p.W_propeller +p.W_ESC + p.W_battery +p.W_wires_power;
% 
% %% Flight Control
% W_flight_control = p.W_receiver + p.n_servos * p.W_servo + p.W_wires_control;
% 
% %External Equipment
% W_external = p.W_landing_gear + p.W_extra_support;

%% Total Weights
W_structure = W_wing + W_HT + W_VT + W_fuselage;
% W_empty = W_structure; %W_powerplant + W_flight_control + W_external ~ if added
% W_takeoff = W_empty;% + p.W_payload ~ if added

%% Weight Breakdown
weight_details.wing.ribs = W_ribs;
weight_details.wing.spars = W_spars;
weight_details.wing.stringers = W_stringers;
weight_details.wing.skin = W_skin;
% weight_details.wing.LE = W_LE;
% weight_details.wing.TE = W_TE;
weight_details.wing.total = W_wing;
    
weight_details.horizontal_tail.total = W_HT;
weight_details.vertical_tail.total = W_VT;
    
weight_details.fuselage.frames = W_frames;
weight_details.fuselage.booms = W_booms;
weight_details.fuselage.longeroons = W_longeroons;
weight_details.fuselage.skin = W_fuse_skin;
weight_details.fuselage.total = W_fuselage;
    
weight_details.structure_total = W_structure;
end