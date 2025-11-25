function [W_wing, W_HT, W_VT, W_fuselage, W_booms, W_structure, W_avionics, W_powerplant, W_total, weight_details, cg_details] = calculate_weights(p)
%% Wing Calculations
AR_wing = p.b_wing^2/p.S_wing;

%Ribs (with lightening factor for cutouts)
rib_chord_avg = p.c_mean;
rib_height = p.tc_wing*rib_chord_avg;
rib_perimeter = 2 *(rib_chord_avg+rib_height); %old
rib_area = rib_perimeter * p.rib_thickness; %old
W_ribs = p.n_ribs * rib_area * p.rib_material;

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
frame_perimeter = 2*(p.W_fuse + p.H_fuse);
W_frames = p.n_frames * frame_perimeter * p.frame_thickness * p.frame_material;

W_longeroons = p.n_longeroons * p.L_fuse * p.longeroon_area * p.longeroon_material;

fuse_surface_area = 2*(p.W_fuse * p.L_fuse + p.H_fuse * p.L_fuse);

W_fuse_skin = fuse_surface_area * p.fuse_skin_material;

W_fuselage = W_frames + W_longeroons + W_fuse_skin;

%% Boom
W_booms = p.n_booms * p.L_boom * p.boom_area * p.boom_material;

%% Powerplant
W_powerplant = 2*p.W_motor + 2*p.W_propeller + 2*p.W_ESC + p.W_battery;%dual prop system

%% Flight Control (Avionics only - excluding powerplant)
W_avionics = p.W_receiver + p.n_servos * p.W_servo;

%% External Equipment
W_external = p.W_landing_gear;

%% Payload
W_payload = p.W_payload;

%% Total Weights
W_structure = W_wing + W_HT + W_VT + W_fuselage + W_booms;
W_empty = W_structure + W_avionics + W_powerplant + W_external;
W_total = W_empty + W_payload;

%% CENTER OF GRAVITY CALCULATIONS FOR EACH PART
% Reference frame: Origin at wing leading edge
% x: positive aft from wing LE
% y: positive right from centerline
% z: positive up from wing LE height

% Wing CG - at 35% chord
cg_details.wing.x_cg = 0.35 * p.c_mean;
cg_details.wing.y_cg = 0;
cg_details.wing.z_cg = 0.5 * p.tc_wing * p.c_mean;
cg_details.wing.weight = W_wing;

% Horizontal Tail CG - at 35% chord
cg_details.HT.x_cg = p.x_HT_LE + 0.35 * c_HT;
cg_details.HT.y_cg = 0;
cg_details.HT.z_cg = p.z_HT + 0.5 * p.HT_plate_thickness;
cg_details.HT.weight = W_HT;

% Vertical Tail CG - at 35% chord, mid-height
cg_details.VT.x_cg = p.x_VT_LE + 0.35 * c_VT;
cg_details.VT.y_cg = 0;
cg_details.VT.z_cg = 0.25 * p.b_VT;
cg_details.VT.weight = W_VT;

% Fuselage CG - at geometric center
cg_details.fuselage.x_cg = p.x_fuse_nose + 0.5 * p.L_fuse;
cg_details.fuselage.y_cg = 0;
cg_details.fuselage.z_cg = p.z_fuse_bottom + 0.5 * p.H_fuse;
cg_details.fuselage.weight = W_fuselage;

% Booms CG - at mid-length
cg_details.booms.x_cg = p.x_boom_start + 0.5 * p.L_fuse;
cg_details.booms.y_cg = 0;
cg_details.booms.z_cg = p.z_boom;
cg_details.booms.weight = W_booms;

%% AVIONICS AND PAYLOAD CG CALCULATIONS

% Motor CG - typically at nose
cg_details.motor.x_cg = p.x_motor;
cg_details.motor.y_cg = 0;
cg_details.motor.z_cg = p.z_motor;
cg_details.motor.weight = 2*p.W_motor; % dual motors

% Propeller CG - forward of motor
cg_details.propeller.x_cg = p.x_propeller;
cg_details.propeller.y_cg = 0;
cg_details.propeller.z_cg = p.z_propeller;
cg_details.propeller.weight = 2*p.W_propeller; % dual props

% ESC CG
cg_details.ESC.x_cg = p.x_ESC;
cg_details.ESC.y_cg = 0;
cg_details.ESC.z_cg = p.z_ESC;
cg_details.ESC.weight = 2*p.W_ESC; % dual ESCs

% Battery CG
cg_details.battery.x_cg = p.x_battery;
cg_details.battery.y_cg = 0;
cg_details.battery.z_cg = p.z_battery;
cg_details.battery.weight = p.W_battery;

% Receiver CG
cg_details.receiver.x_cg = p.x_receiver;
cg_details.receiver.y_cg = 0;
cg_details.receiver.z_cg = p.z_receiver;
cg_details.receiver.weight = p.W_receiver;

% Servos CG 
cg_details.servos.x_cg = p.x_servos;
cg_details.servos.y_cg = 0;
cg_details.servos.z_cg = p.z_servos;
cg_details.servos.weight = p.n_servos * p.W_servo;

% Landing Gear CG
cg_details.landing_gear.x_cg = p.x_landing_gear;
cg_details.landing_gear.y_cg = 0;
cg_details.landing_gear.z_cg = p.z_landing_gear;
cg_details.landing_gear.weight = p.W_landing_gear;

% Payload CG
cg_details.payload.x_cg = p.x_payload;
cg_details.payload.y_cg = 0;
cg_details.payload.z_cg = p.z_payload;
cg_details.payload.weight = W_payload;

%% Weight Breakdown
weight_details.wing.ribs = W_ribs;
weight_details.wing.spars = W_spars;
weight_details.wing.stringers = W_stringers;
weight_details.wing.skin = W_skin;
weight_details.wing.total = W_wing;
    
weight_details.horizontal_tail.total = W_HT;
weight_details.vertical_tail.total = W_VT;
    
weight_details.fuselage.frames = W_frames;
weight_details.fuselage.booms = W_booms;
weight_details.fuselage.longeroons = W_longeroons;
weight_details.fuselage.skin = W_fuse_skin;
weight_details.fuselage.total = W_fuselage;

weight_details.powerplant.motor = 2*p.W_motor;
weight_details.powerplant.propeller = 2*p.W_propeller;
weight_details.powerplant.ESC = 2*p.W_ESC;
weight_details.powerplant.battery = p.W_battery;
weight_details.powerplant.total = W_powerplant;

weight_details.flight_control.receiver = p.W_receiver;
weight_details.flight_control.servos = p.n_servos * p.W_servo;
weight_details.flight_control.total = W_avionics;

weight_details.external.landing_gear = p.W_landing_gear;
weight_details.external.extra_support = p.W_extra_support;
weight_details.external.total = W_external;

weight_details.payload = W_payload;
    
weight_details.structure_total = W_structure;
weight_details.avionics_total = W_avionics;
weight_details.powerplant_total = W_powerplant;
weight_details.empty_weight = W_empty;
weight_details.takeoff_weight = W_total;

%% TOTAL AIRCRAFT CG CALCULATION (Empty Weight)
W_empty_total = W_structure + W_avionics + W_powerplant + W_external;

x_cg_empty = (cg_details.wing.weight * cg_details.wing.x_cg + ...
              cg_details.HT.weight * cg_details.HT.x_cg + ...
              cg_details.VT.weight * cg_details.VT.x_cg + ...
              cg_details.fuselage.weight * cg_details.fuselage.x_cg + ...
              cg_details.booms.weight * cg_details.booms.x_cg + ...
              cg_details.motor.weight * cg_details.motor.x_cg + ...
              cg_details.propeller.weight * cg_details.propeller.x_cg + ...
              cg_details.ESC.weight * cg_details.ESC.x_cg + ...
              cg_details.battery.weight * cg_details.battery.x_cg + ...
              cg_details.receiver.weight * cg_details.receiver.x_cg + ...
              cg_details.servos.weight * cg_details.servos.x_cg + ...
              cg_details.landing_gear.weight * cg_details.landing_gear.x_cg) / W_empty_total;

y_cg_empty = 0; % Should be 0 for symmetric aircraft

z_cg_empty = (cg_details.wing.weight * cg_details.wing.z_cg + ...
              cg_details.HT.weight * cg_details.HT.z_cg + ...
              cg_details.VT.weight * cg_details.VT.z_cg + ...
              cg_details.fuselage.weight * cg_details.fuselage.z_cg + ...
              cg_details.booms.weight * cg_details.booms.z_cg + ...
              cg_details.motor.weight * cg_details.motor.z_cg + ...
              cg_details.propeller.weight * cg_details.propeller.z_cg + ...
              cg_details.ESC.weight * cg_details.ESC.z_cg + ...
              cg_details.battery.weight * cg_details.battery.z_cg + ...
              cg_details.receiver.weight * cg_details.receiver.z_cg + ...
              cg_details.servos.weight * cg_details.servos.z_cg + ...
              cg_details.landing_gear.weight * cg_details.landing_gear.z_cg) / W_empty_total;

% Store empty weight CG
cg_details.aircraft_empty.x_cg = x_cg_empty;
cg_details.aircraft_empty.y_cg = y_cg_empty;
cg_details.aircraft_empty.z_cg = z_cg_empty;
cg_details.aircraft_empty.weight = W_empty_total;

%% TOTAL AIRCRAFT CG CALCULATION (Takeoff Weight with Payload)
x_cg_total = (W_empty_total * x_cg_empty + ...
              cg_details.payload.weight * cg_details.payload.x_cg) / W_total;

z_cg_total = (W_empty_total * z_cg_empty + ...
              cg_details.payload.weight * cg_details.payload.z_cg) / W_total;

% Store takeoff weight CG
cg_details.aircraft_takeoff.x_cg = x_cg_total;
cg_details.aircraft_takeoff.y_cg = 0;
cg_details.aircraft_takeoff.z_cg = z_cg_total;
cg_details.aircraft_takeoff.weight = W_total;

end