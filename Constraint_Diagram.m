% Assumed
CLmax = 1.0;
CL_TO = CLmax * 0.8;
e = 0.7;
AR = 3.6;
CD0 = 0.12;


% Speeds (ft/s)
Vinf   = 27.5 * 3.28084    % 0.3048
Vstall = Vinf * (2/3)      % 2/3  
V_TO   = 1.2 * Vstall   


% Distance
Sg = 150;       


% Rolling friction coefficient
mu = 0.02;
CD_TO = 0.2;


% ROC
d1 = 500 - (Vinf * 1);
t1 = d1 / V_TO;
ROC = 100 / t1;
% ROC_fpm = 500;           % ft/min
% ROC = ROC_fpm * 0.00508; % 2.54 m/s


% Load factor
n = 2;


% Constants
g = 32.1740;                    % ft/s^2
rho = 0.0023062065;             % slugs/ft^3
MTOW = 6.2;

% Induced drag factor
K = 1 / (pi * e * AR);


% Dynamic pressure
q_cruise = 0.5 * rho * Vinf^2;
q_to     = 0.5 * rho * V_TO^2;

%Thrust
CDi = 0.35^2 * K;
T = 0.5 * 1.198 * 27.5^2 * 0.2322576 * (CD0 + CDi) * 1.25

% W/S sweep (N/m^2)
WS = linspace(0.001, 8, 50000);


% Stall constraint 
WS_stall = 0.5 * rho * Vstall^2 * CLmax;
Sreal = MTOW / 3.5
b = 3;
ARreal = b^2 / Sreal
AR2 = 5;
Sreal2 = 1/(ARreal/(b^2))

% Takeoff constraint
TW_takeoff = (V_TO^2) / (2*g*Sg) ...
         + (q_to * CD_TO) ./ WS ...
         + mu * (1 - (q_to * CL_TO) ./ WS);


% Cruise constraint
TW_cruise = (q_cruise * CD0) ./ WS ...
        + WS ./ (q_cruise * pi * e * AR);


% ROC constraint
TW_ROC = ROC / Vinf ...
     + (q_cruise * CD0) ./ WS ...
     + K * WS ./ q_cruise;


% Maximum load factor
TW_turn = q_cruise * ((CD0./WS) + (K .* ( (n.^2 ./ q_cruise.^2) ) .* WS));


% Plot
color_pass  = [0, 0, 0] / 255;   
color_laps  = [172, 43, 55] / 255;   
color_cargo = [169, 176, 183] / 255; 
color_EF    = [195, 17, 47] / 255;

figure; hold on; grid on;
plot(WS, TW_takeoff, 'Color', color_pass, 'LineWidth', 1, 'DisplayName','Takeoff');
plot(WS, 1.3*TW_cruise, 'Color', color_laps, 'MarkerEdgeColor', color_laps, 'LineWidth', 1, 'DisplayName','Cruise');
plot(WS, TW_ROC, 'Color', color_cargo, 'LineWidth', 1, 'DisplayName','ROC');
plot(WS, 1.3*TW_turn, 'Color', color_EF, 'LineWidth', 1, 'DisplayName','Turn constraint');
xline(WS_stall, '--', 'Color', color_pass, 'DisplayName','Stall Constraint', 'LineWidth',2);
xlabel('W/S (lb/ft^2)');
ylabel('T/W');
title('Constraint Diagram with Design Point');
legend('Takeoff', 'Cruise', 'ROC', 'Turn', 'Stall')
ylim([0 3.5]);
xlim([0 8]);

% add historical planes
tw_mich = 0.94;
ws_mich = 3.57;
tw_next = 0.726;
ws_next = 3.58;
plot(ws_mich, tw_mich, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Michigan (2020)');
%text(ws_mich-0.5, tw_mich + 0.1, 'Michigan (2021)', 'FontSize',6);
plot(ws_next, tw_next, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Dayananda (2020)');
%text(ws_next-0.5, tw_next + 0.1, 'Dayananda (2021)', 'FontSize',6);
plot(3.16, 1.425, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','VT (2024)');
%%text(3.16, 1.425 + 0.1, 'VT (2024)', 'FontSize',6);
plot(2.603, 2.108590554, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','LV (2020)');
%text(2.603+0.1, 2.108590554 + 0.1, 'LV (2020)', 'FontSize',6);
plot(0.98,1.529411765, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','GIT (2020)');
%text(0.98-0.5, 1.529411765 + 0.1, 'GIT (2020)', 'FontSize',6);
plot(1.38,0.641, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','USC (2020)');
%text(1.38-0.5, 0.641 + 0.1, 'USC (2020)', 'FontSize',6);
plot(5.143,0.4, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UNSW (2025)');
%text(5.143-0.5,0.4 + 0.1, 'UNSW (2025)', 'FontSize',6);
plot(1.733,2.113, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','USC (2024)');
%%text(1.733,2.113 + 0.1, 'USC (2024)', 'FontSize',6);
plot(1.487,1.6, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','WU St Louis (2022)');
%%text(1.487,1.6 + 0.1, 'WU St Louis (2022)', 'FontSize',6);
plot(0.953,0.889, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','In Thrust We Trust (2018)');
%text(0.953-0.5,0.889 + 0.1, 'In Thrust We Trust (2018)', 'FontSize',6);
plot(0.965,0.829, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Sequoia (2025)');
%text(0.965-0.5,0.829 + 0.1, 'Sequoia (2025)', 'FontSize',6);
plot(0.846,0.912, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Tri Spectra (2018)');
%text(0.846-0.5,0.912 + 0.1, 'Tri Spectra (2018)', 'Fontsize',6);
plot(3.722,0.481, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','El Agave (2019)');
%text(3.722-0.5,0.481 + 0.1, 'El Agave (2019)', 'FontSize',6);
plot(0.338,0.686, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Freddy (2022)');
%text(0.338-0.5,0.686 + 0.1, 'Freddy (2022)', 'FontSize',6);
plot(0.5072,0.599, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','REVADA (2022)');
%text(0.5072-0.5,0.599 + 0.1, 'Revada (2022)', 'FontSize',6);
plot(1.167,0.699, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Robot Chicken (2024)');
%text(1.167-0.5,0.699 + 0.1, 'Robot Chicken (2024)', 'FontSize',6);
%plot(1.032,0.963, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','AEMS STOL (2024)');
%text(1.032-0.5,0.963 + 0.1, 'AEMS STOL (2024)', 'FontSize',6);
%plot(1.016,1.508, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','AeroLink STOL (2024)');
%text(1.016+0.1,1.508 - 0.05, 'AeroLink STOL (2024)', 'FontSize',6);
%plot(0.945,1.3, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Deuces STOL (2024)');
%text(0.945+0.2,1.3 -0.05, 'Deuces STOL (2024)', 'FontSize',6);
%plot(1,0.865, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Starfleet (2024)');
%text(1-0.5,0.865 + 0.1, 'Starfleet STOL (2024)', 'FontSize',6);
%plot(0.465,0.523, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Sky-King');
%text(0.465-0.5,0.523 + 0.1, 'Sky-King', 'FontSize',6);
%plot(1.615,0.45, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Marlin Arrow');
%text(1.615-0.5,0.45 + 0.1, 'Marlin Arrow', 'FontSize',6);
%plot(0.251,0.819, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Hobby Zone');
%text(0.251-0.5,0.819 + 0.1, 'Hobby Zone', 'FontSize',6);
%plot(0.308,0.313, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','WL F-949');
%text(0.308-0.5,0.313 + 0.1, 'WL F-949', 'FontSize',6);
%plot(0.5,0.367, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','XK A-800');
%text(0.5-0.5,0.367 + 0.1, 'XK A-800', 'FontSize',6);
plot(3.9,0.845, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UCLA (2025)');
%text(3.9-0.5,0.845 + 0.1, 'UCLA (2025)', 'FontSize',6);
plot(3.5, 1, 'kp', 'MarkerSize',12,'MarkerFaceColor','y', 'DisplayName','UCLA (2025)');
text(4.25, 1, 'Design Point', 'FontSize',10, 'FontWeight','bold');
legend('Takeoff', 'Cruise', 'ROC', 'Turn', 'Stall');


%plot(182.6609,1.2465, 'y*', 'MarkerSize',10,'Color',, color_laps, 'HandleVisibility','off');
% saveas(gcf, fullfile('C:\Users\gmbol\OneDrive\School Work\MQP', 'Constraint Diagram Design Pt.png'));

% Sensitivity 1

figure; hold on; grid on;

% Baseline plots
plot(WS, TW_takeoff, 'b-', 'LineWidth', 1.5, 'DisplayName','Takeoff');
plot(WS, 1.3*TW_cruise,  'm-', 'LineWidth', 1.5, 'DisplayName','Cruise');
plot(WS, 1.3*TW_ROC,     'g-', 'LineWidth', 1.5, 'DisplayName','ROC');
plot(WS, 1.3*TW_turn,    'c-', 'LineWidth', 1.5, 'DisplayName','Turn');
xline(WS_stall, '--r', 'LineWidth',1.5, 'DisplayName','Stall');

% Sensitivity variables
CLmax_vals = [0.8*CLmax, 1.2*CLmax];
Sg_vals    = [0.8*Sg,    1.2*Sg];
zeta_vals  = [0.8, 1.2];
ROC_vals   = [0.8*ROC,   1.2*ROC];
n_vals     = [0.8*n,     1.2*n];

% --- Stall (CLmax variation) ---
for CLm = CLmax_vals
    WS_stall_sens = 0.5 * rho * Vstall^2 * CLm;
    xline(WS_stall_sens, 'r--', 'HandleVisibility','off');
end

% --- Takeoff (Sg variation) ---
for Sg_ = Sg_vals
    %WS = 182.6609 N/M^2;
    TW_takeoff_sens = (V_TO^2) / (2*g*Sg_) ...
         + (q_to * CD_TO) ./ WS ...
         + mu * (1 - (q_to * CL_TO) ./ WS);
    plot(WS, TW_takeoff_sens, 'b--', 'HandleVisibility','off');
end

% --- Cruise (ζ variation) ---
for zeta = zeta_vals
    TW_cruise_sens = (q_cruise * CD0) ./ WS ...
        + zeta * WS ./ (q_cruise * pi * e * AR);
    plot(WS, 1.3*TW_cruise_sens, 'm--', 'HandleVisibility','off');
end

% --- ROC (ROC variation) ---
for ROC_ = ROC_vals
    TW_ROC_sens = ROC_ / Vinf ...
     + (q_cruise * CD0) ./ WS ...
     + K * WS ./ q_cruise;
    plot(WS, 1.3*TW_ROC_sens, 'g--', 'HandleVisibility','off');
end

% --- Turn (n variation) ---
for n_ = n_vals
    TW_turn_sens = (q_cruise * CD0) ./ WS ...
                  + K * n_^2 .* WS ./ q_cruise;
    plot(WS, 1.3*TW_turn_sens, 'c--', 'HandleVisibility','off');
end

xlabel('W/S (lb/ft^2)');
ylabel('T/W');
title('Constraint Diagram with Sensitivity (20%)');
ylim([0 3.5]);
xlim([0 8]);

% add historical planes
tw_mich = 0.94;
ws_mich = 3.57;
tw_next = 0.726;
ws_next = 3.58;
plot(ws_mich, tw_mich, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Michigan (2020)');
%text(ws_mich-0.5, tw_mich + 0.1, 'Michigan (2021)', 'FontSize',6);
plot(ws_next, tw_next, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Dayananda (2020)');
%text(ws_next-0.5, tw_next + 0.1, 'Dayananda (2021)', 'FontSize',6);
plot(3.16, 1.425, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','VT (2024)');
%%text(3.16, 1.425 + 0.1, 'VT (2024)', 'FontSize',6);
plot(2.603, 2.108590554, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','LV (2020)');
%text(2.603+0.1, 2.108590554 + 0.1, 'LV (2020)', 'FontSize',6);
plot(0.98,1.529411765, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','GIT (2020)');
%text(0.98-0.5, 1.529411765 + 0.1, 'GIT (2020)', 'FontSize',6);
plot(1.38,0.641, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','USC (2020)');
%text(1.38-0.5, 0.641 + 0.1, 'USC (2020)', 'FontSize',6);
plot(5.143,0.4, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UNSW (2025)');
%text(5.143-0.5,0.4 + 0.1, 'UNSW (2025)', 'FontSize',6);
plot(1.733,2.113, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','USC (2024)');
%%text(1.733,2.113 + 0.1, 'USC (2024)', 'FontSize',6);
plot(1.487,1.6, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','WU St Louis (2022)');
%%text(1.487,1.6 + 0.1, 'WU St Louis (2022)', 'FontSize',6);
plot(0.953,0.889, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','In Thrust We Trust (2018)');
%text(0.953-0.5,0.889 + 0.1, 'In Thrust We Trust (2018)', 'FontSize',6);
plot(0.965,0.829, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Sequoia (2025)');
%text(0.965-0.5,0.829 + 0.1, 'Sequoia (2025)', 'FontSize',6);
plot(0.846,0.912, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Tri Spectra (2018)');
%text(0.846-0.5,0.912 + 0.1, 'Tri Spectra (2018)', 'Fontsize',6);
plot(3.722,0.481, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','El Agave (2019)');
%text(3.722-0.5,0.481 + 0.1, 'El Agave (2019)', 'FontSize',6);
plot(0.338,0.686, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Freddy (2022)');
%text(0.338-0.5,0.686 + 0.1, 'Freddy (2022)', 'FontSize',6);
plot(0.5072,0.599, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','REVADA (2022)');
%text(0.5072-0.5,0.599 + 0.1, 'Revada (2022)', 'FontSize',6);
plot(1.167,0.699, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Robot Chicken (2024)');
%text(1.167-0.5,0.699 + 0.1, 'Robot Chicken (2024)', 'FontSize',6);
%plot(1.032,0.963, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','AEMS STOL (2024)');
%text(1.032-0.5,0.963 + 0.1, 'AEMS STOL (2024)', 'FontSize',6);
%plot(1.016,1.508, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','AeroLink STOL (2024)');
%text(1.016+0.1,1.508 - 0.05, 'AeroLink STOL (2024)', 'FontSize',6);
%plot(0.945,1.3, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Deuces STOL (2024)');
%text(0.945+0.2,1.3 -0.05, 'Deuces STOL (2024)', 'FontSize',6);
%plot(1,0.865, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Starfleet (2024)');
%text(1-0.5,0.865 + 0.1, 'Starfleet STOL (2024)', 'FontSize',6);
%plot(0.465,0.523, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Sky-King');
%text(0.465-0.5,0.523 + 0.1, 'Sky-King', 'FontSize',6);
%plot(1.615,0.45, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Marlin Arrow');
%text(1.615-0.5,0.45 + 0.1, 'Marlin Arrow', 'FontSize',6);
%plot(0.251,0.819, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Hobby Zone');
%text(0.251-0.5,0.819 + 0.1, 'Hobby Zone', 'FontSize',6);
%plot(0.308,0.313, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','WL F-949');
%text(0.308-0.5,0.313 + 0.1, 'WL F-949', 'FontSize',6);
%plot(0.5,0.367, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','XK A-800');
%text(0.5-0.5,0.367 + 0.1, 'XK A-800', 'FontSize',6);
plot(3.9,0.845, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UCLA (2025)');
%text(3.9-0.5,0.845 + 0.1, 'UCLA (2025)', 'FontSize',6);
plot(4, 1, 'kp', 'MarkerSize',12,'Color','y', 'DisplayName','UCLA (2025)');
text(4.25, 1, 'Design Point', 'FontSize',10, 'FontWeight','bold');
legend('Takeoff', 'Cruise', 'ROC', 'Turn', 'Stall');



%plot(182.6609,1.2465, 'y*', 'MarkerSize',10,'Color',, color_laps, 'HandleVisibility','off');
saveas(gcf, fullfile('C:\Users\gmbol\OneDrive\School Work\MQP', 'Constraint Diagram First Sensitivity.png'));

%% ---- Updated Helper function for constraints ----
compute_constraints = @(Vstall_val, CD0_val) deal( ...
    0.5 * rho * Vstall_val^2 * CLmax, ...                           % WS_stall
    ( (1.2*Vstall_val)^2 )/(2*g*Sg) + (q_to*CD_TO)./WS + mu*(1 - (q_to*CL_TO)./WS), ... % Takeoff (V_TO = 1.2*Vstall)
    (0.5*rho*(1.5*Vstall_val)^2*CD0_val)./WS + WS./(0.5*rho*(1.5*Vstall_val)^2*pi*e*AR), ... % Cruise (Vcruise = 1.5*Vstall)
    ROC/(1.5*Vstall_val) + (0.5*rho*(1.5*Vstall_val)^2*CD0_val)./WS + K*WS./(0.5*rho*(1.5*Vstall_val)^2), ... % ROC
    (0.5*rho*(1.5*Vstall_val)^2*CD0_val)./WS + K*n^2.*WS./(0.5*rho*(1.5*Vstall_val)^2) );   % Turn

%% ---- Baseline constraints ----
[WS_stall, TW_takeoff, TW_cruise, TW_ROC, TW_turn] = compute_constraints(Vstall, CD0);

%% ---- Sensitivity (20%) ----
Vstall_low  = 0.8 * Vstall;
Vstall_high = 1.2 * Vstall;
CD0_low     = 0.8 * CD0;
CD0_high    = 1.2 * CD0;

[WS_stall_low,  TO_low,  CR_low,  ROC_low,  TURN_low]  = compute_constraints(Vstall_low,  CD0);
[WS_stall_high, TO_high, CR_high, ROC_high, TURN_high] = compute_constraints(Vstall_high, CD0);
[~,             TO_CDlow,  CR_CDlow,  ROC_CDlow,  TURN_CDlow]  = compute_constraints(Vstall, CD0_low);
[~,             TO_CDhigh, CR_CDhigh, ROC_CDhigh, TURN_CDhigh] = compute_constraints(Vstall, CD0_high);

%% ---- Plot 1: Baseline + Vstall sensitivity ----
figure; hold on; grid on;

% Baseline
plot(WS, TW_takeoff, 'b-', 'LineWidth', 1.5, 'DisplayName','Takeoff (baseline)');
plot(WS, 1.3*TW_cruise,  'm-', 'LineWidth', 1.5, 'DisplayName','Cruise (baseline)');
plot(WS, 1.3*TW_ROC,     'g-', 'LineWidth', 1.5, 'DisplayName','ROC (baseline)');
plot(WS, 1.3*TW_turn,    'c-', 'LineWidth', 1.5, 'DisplayName','Turn (baseline)');
xline(WS_stall, 'r--', 'LineWidth', 1.5, 'DisplayName','Stall (baseline)');

% Vstall sensitivity
plot(WS, TO_low,  'b--', 'DisplayName','Takeoff (Vstall -20%)');
plot(WS, TO_high, 'b--', 'HandleVisibility','off');
plot(WS, 1.3*CR_low,  'm--', 'DisplayName','Cruise (Vstall -20%)');
plot(WS, 1.3*CR_high, 'm--', 'HandleVisibility','off');
plot(WS, 1.3*ROC_low, 'g--', 'DisplayName','ROC (Vstall -20%)');
plot(WS, 1.3*ROC_high,'g--', 'HandleVisibility','off');
plot(WS, 1.3*TURN_low,'c--', 'DisplayName','Turn (Vstall -20%)');
plot(WS, 1.3*TURN_high,'c--','HandleVisibility','off');
xline(WS_stall_low,  'r--', 'DisplayName','Stall (Vstall -20%)');
xline(WS_stall_high, 'r--', 'HandleVisibility','off');

% add historical planes
tw_mich = 0.94;
ws_mich = 3.57;
tw_next = 0.726;
ws_next = 3.58;
plot(ws_mich, tw_mich, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Michigan (2020)');
%text(ws_mich-0.5, tw_mich + 0.1, 'Michigan (2021)', 'FontSize',6);
plot(ws_next, tw_next, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Dayananda (2020)');
%text(ws_next-0.5, tw_next + 0.1, 'Dayananda (2021)', 'FontSize',6);
plot(3.16, 1.425, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','VT (2024)');
%%text(3.16, 1.425 + 0.1, 'VT (2024)', 'FontSize',6);
plot(2.603, 2.108590554, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','LV (2020)');
%text(2.603+0.1, 2.108590554 + 0.1, 'LV (2020)', 'FontSize',6);
plot(0.98,1.529411765, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','GIT (2020)');
%text(0.98-0.5, 1.529411765 + 0.1, 'GIT (2020)', 'FontSize',6);
plot(1.38,0.641, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','USC (2020)');
%text(1.38-0.5, 0.641 + 0.1, 'USC (2020)', 'FontSize',6);
plot(5.143,0.4, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UNSW (2025)');
%text(5.143-0.5,0.4 + 0.1, 'UNSW (2025)', 'FontSize',6);
plot(1.733,2.113, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','USC (2024)');
%%text(1.733,2.113 + 0.1, 'USC (2024)', 'FontSize',6);
plot(1.487,1.6, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','WU St Louis (2022)');
%%text(1.487,1.6 + 0.1, 'WU St Louis (2022)', 'FontSize',6);
plot(0.953,0.889, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','In Thrust We Trust (2018)');
%text(0.953-0.5,0.889 + 0.1, 'In Thrust We Trust (2018)', 'FontSize',6);
plot(0.965,0.829, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Sequoia (2025)');
%text(0.965-0.5,0.829 + 0.1, 'Sequoia (2025)', 'FontSize',6);
plot(0.846,0.912, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Tri Spectra (2018)');
%text(0.846-0.5,0.912 + 0.1, 'Tri Spectra (2018)', 'Fontsize',6);
plot(3.722,0.481, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','El Agave (2019)');
%text(3.722-0.5,0.481 + 0.1, 'El Agave (2019)', 'FontSize',6);
plot(0.338,0.686, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Freddy (2022)');
%text(0.338-0.5,0.686 + 0.1, 'Freddy (2022)', 'FontSize',6);
plot(0.5072,0.599, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','REVADA (2022)');
%text(0.5072-0.5,0.599 + 0.1, 'Revada (2022)', 'FontSize',6);
plot(1.167,0.699, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Robot Chicken (2024)');
%text(1.167-0.5,0.699 + 0.1, 'Robot Chicken (2024)', 'FontSize',6);
%plot(1.032,0.963, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','AEMS STOL (2024)');
%text(1.032-0.5,0.963 + 0.1, 'AEMS STOL (2024)', 'FontSize',6);
%plot(1.016,1.508, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','AeroLink STOL (2024)');
%text(1.016+0.1,1.508 - 0.05, 'AeroLink STOL (2024)', 'FontSize',6);
%plot(0.945,1.3, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Deuces STOL (2024)');
%text(0.945+0.2,1.3 -0.05, 'Deuces STOL (2024)', 'FontSize',6);
%plot(1,0.865, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Starfleet (2024)');
%text(1-0.5,0.865 + 0.1, 'Starfleet STOL (2024)', 'FontSize',6);
%plot(0.465,0.523, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Sky-King');
%text(0.465-0.5,0.523 + 0.1, 'Sky-King', 'FontSize',6);
%plot(1.615,0.45, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Marlin Arrow');
%text(1.615-0.5,0.45 + 0.1, 'Marlin Arrow', 'FontSize',6);
%plot(0.251,0.819, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Hobby Zone');
%text(0.251-0.5,0.819 + 0.1, 'Hobby Zone', 'FontSize',6);
%plot(0.308,0.313, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','WL F-949');
%text(0.308-0.5,0.313 + 0.1, 'WL F-949', 'FontSize',6);
%plot(0.5,0.367, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','XK A-800');
%text(0.5-0.5,0.367 + 0.1, 'XK A-800', 'FontSize',6);
plot(3.9,0.845, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UCLA (2025)');
%text(3.9-0.5,0.845 + 0.1, 'UCLA (2025)', 'FontSize',6);
plot(4, 1, 'kp', 'MarkerSize',12,'Color','y', 'DisplayName','UCLA (2025)');
text(4.25, 1, 'Design Point', 'FontSize',10, 'FontWeight','bold');
legend('Takeoff', 'Cruise', 'ROC', 'Turn', 'Stall');


xlabel('W/S (lb/ft^2)');
ylabel('T/W');
title('Constraint Diagram – V_{stall} Sensitivity (±20%)');
legend('Location','northeast');
ylim([0 3.5]); xlim([0 8]);
saveas(gcf, fullfile('C:\Users\gmbol\OneDrive\School Work\MQP', 'Constraint Diagram Vstall Sensitivity.png'));


%% ---- Plot 2: Baseline + CD0 sensitivity ----
figure; hold on; grid on;

% Baseline
plot(WS, TW_takeoff, 'b-', 'LineWidth', 1.5, 'DisplayName','Takeoff (baseline)');
plot(WS, 1.3*TW_cruise,  'm-', 'LineWidth', 1.5, 'DisplayName','Cruise (baseline)');
plot(WS, 1.3*TW_ROC,     'g-', 'LineWidth', 1.5, 'DisplayName','ROC (baseline)');
plot(WS, 1.3*TW_turn,    'c-', 'LineWidth', 1.5, 'DisplayName','Turn (baseline)');
xline(WS_stall, 'r--', 'LineWidth', 1.5, 'DisplayName','Stall (baseline)');

% CD0 sensitivity
plot(WS, TO_CDlow,  'b:', 'DisplayName','Takeoff (CD0 -20%)');
plot(WS, TO_CDhigh, 'b:', 'HandleVisibility','off');
plot(WS, 1.3*CR_CDlow,  'm:', 'DisplayName','Cruise (CD0 -20%)');
plot(WS, 1.3*CR_CDhigh, 'm:', 'HandleVisibility','off');
plot(WS, 1.3*ROC_CDlow, 'g:', 'DisplayName','ROC (CD0 -20%)');
plot(WS, 1.3*ROC_CDhigh,'g:', 'HandleVisibility','off');
plot(WS, 1.3*TURN_CDlow,'c:', 'DisplayName','Turn (CD0 -20%)');
plot(WS, 1.3*TURN_CDhigh,'c:','HandleVisibility','off');

% add historical planes
tw_mich = 0.94;
ws_mich = 3.57;
tw_next = 0.726;
ws_next = 3.58;
plot(ws_mich, tw_mich, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Michigan (2020)');
%text(ws_mich-0.5, tw_mich + 0.1, 'Michigan (2021)', 'FontSize',6);
plot(ws_next, tw_next, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Dayananda (2020)');
%text(ws_next-0.5, tw_next + 0.1, 'Dayananda (2021)', 'FontSize',6);
plot(3.16, 1.425, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','VT (2024)');
%%text(3.16, 1.425 + 0.1, 'VT (2024)', 'FontSize',6);
plot(2.603, 2.108590554, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','LV (2020)');
%text(2.603+0.1, 2.108590554 + 0.1, 'LV (2020)', 'FontSize',6);
plot(0.98,1.529411765, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','GIT (2020)');
%text(0.98-0.5, 1.529411765 + 0.1, 'GIT (2020)', 'FontSize',6);
plot(1.38,0.641, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','USC (2020)');
%text(1.38-0.5, 0.641 + 0.1, 'USC (2020)', 'FontSize',6);
plot(5.143,0.4, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UNSW (2025)');
%text(5.143-0.5,0.4 + 0.1, 'UNSW (2025)', 'FontSize',6);
plot(1.733,2.113, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','USC (2024)');
%%text(1.733,2.113 + 0.1, 'USC (2024)', 'FontSize',6);
plot(1.487,1.6, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','WU St Louis (2022)');
%%text(1.487,1.6 + 0.1, 'WU St Louis (2022)', 'FontSize',6);
plot(0.953,0.889, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','In Thrust We Trust (2018)');
%text(0.953-0.5,0.889 + 0.1, 'In Thrust We Trust (2018)', 'FontSize',6);
plot(0.965,0.829, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Sequoia (2025)');
%text(0.965-0.5,0.829 + 0.1, 'Sequoia (2025)', 'FontSize',6);
plot(0.846,0.912, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Tri Spectra (2018)');
%text(0.846-0.5,0.912 + 0.1, 'Tri Spectra (2018)', 'Fontsize',6);
plot(3.722,0.481, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','El Agave (2019)');
%text(3.722-0.5,0.481 + 0.1, 'El Agave (2019)', 'FontSize',6);
plot(0.338,0.686, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Freddy (2022)');
%text(0.338-0.5,0.686 + 0.1, 'Freddy (2022)', 'FontSize',6);
plot(0.5072,0.599, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','REVADA (2022)');
%text(0.5072-0.5,0.599 + 0.1, 'Revada (2022)', 'FontSize',6);
plot(1.167,0.699, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','Robot Chicken (2024)');
%text(1.167-0.5,0.699 + 0.1, 'Robot Chicken (2024)', 'FontSize',6);
%plot(1.032,0.963, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','AEMS STOL (2024)');
%text(1.032-0.5,0.963 + 0.1, 'AEMS STOL (2024)', 'FontSize',6);
%plot(1.016,1.508, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','AeroLink STOL (2024)');
%text(1.016+0.1,1.508 - 0.05, 'AeroLink STOL (2024)', 'FontSize',6);
%plot(0.945,1.3, 'gp', 'MarkerSize',6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName','Deuces STOL (2024)');
%text(0.945+0.2,1.3 -0.05, 'Deuces STOL (2024)', 'FontSize',6);
%plot(1,0.865, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Starfleet (2024)');
%text(1-0.5,0.865 + 0.1, 'Starfleet STOL (2024)', 'FontSize',6);
%plot(0.465,0.523, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Sky-King');
%text(0.465-0.5,0.523 + 0.1, 'Sky-King', 'FontSize',6);
%plot(1.615,0.45, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Marlin Arrow');
%text(1.615-0.5,0.45 + 0.1, 'Marlin Arrow', 'FontSize',6);
%plot(0.251,0.819, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','Hobby Zone');
%text(0.251-0.5,0.819 + 0.1, 'Hobby Zone', 'FontSize',6);
%plot(0.308,0.313, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','WL F-949');
%text(0.308-0.5,0.313 + 0.1, 'WL F-949', 'FontSize',6);
%plot(0.5,0.367, 'p', 'MarkerSize',6,'Color',, color_laps, 'DisplayName','XK A-800');
%text(0.5-0.5,0.367 + 0.1, 'XK A-800', 'FontSize',6);
plot(3.9,0.845, 'p', 'MarkerSize',6,'MarkerFaceColor', color_laps, 'MarkerEdgeColor', color_laps, 'DisplayName','UCLA (2025)');
%text(3.9-0.5,0.845 + 0.1, 'UCLA (2025)', 'FontSize',6);
plot(4, 1, 'kp', 'MarkerSize',12,'Color','y', 'DisplayName','UCLA (2025)');
text(4.25, 1, 'Design Point', 'FontSize',10, 'FontWeight','bold');
legend('Takeoff', 'Cruise', 'ROC', 'Turn', 'Stall');


xlabel('W/S (lb/ft^2)');
ylabel('T/W');
title('Constraint Diagram – C_{D0} Sensitivity (±20%)');
legend('Location','northeast');
ylim([0 3.5]); xlim([0 8]);
saveas(gcf, fullfile('C:\Users\gmbol\OneDrive\School Work\MQP', 'Constraint Diagram CD0.png'));



