% V-n

g = 9.81;
c = 10/12;
WS = 3.5;                 % lb/ft^2
rho = 1.198 * 0.00194032; % slugs/ft^3
Cl_max = 1.0;
Cn_max = Cl_max*1.1;
CL3 = 1.6;
Kc = 33;                  
Ude = 12;
Ude2 = 6;                 % fps
V = linspace(0,250, 100); 

% Kg and mug
mug = 2 * (WS / (rho * c * g * CL3));
Kg = 0.88 * (mug / (5.3 + mug));

% nlims
nlimpos = 3;
nlimneg = -0.4 * nlimpos;

% Stall and Cruise
Vstall = sqrt(2 * (WS / (rho * Cn_max)));
Vcruise_min = Kc*sqrt(WS);
Vdive_min = 1.25 * Vcruise_min;
Va = Vstall * sqrt(nlimpos);

% Stall lines
n = 0.5 .* rho .* V.^2 .* Cn_max .* (1/WS);
n_neg = -n;

% Gust lines
nlimgustC = 1 + ((Kg.*Ude.*V.*CL3)./(498 .* WS));
nlimgustCneg = 1 - ((Kg.*Ude.*V.*CL3)./(498 .* WS));
nlimgustD = 1 + ((Kg.*Ude2.*V.*CL3)./(498 .* WS));
nlimgustDneg = 1 - ((Kg.*Ude2.*V.*CL3)./(498 .* WS));

% Plot 
figure; hold on; grid on;
plot(V, n, 'b', 'LineWidth', 1);
plot(V, n_neg, 'b','LineWidth', 1);
yline(nlimpos, '--', 'LineWidth', 1);
yline(nlimneg, '--', 'LineWidth', 1);
plot(V, nlimgustC, '--r', 'LineWidth', 1);
plot(V, nlimgustCneg, '--r', 'LineWidth', 1);
plot(V, nlimgustD, '--g', 'LineWidth', 1);
plot(V, nlimgustDneg, '--g', 'LineWidth', 1);
plot(Vcruise_min, 0, 'yo', 'MarkerFaceColor','y');
plot(Vstall, 0, 'ko', 'MarkerFaceColor', 'k');
plot(Vdive_min, 0, 'bo', 'MarkerFaceColor', 'b');
%plot(Va, 0, 'ko', 'MarkerFaceColor', 'k');
title('V-n Diagram');
xlabel('V (ft/s)');
ylabel('n (Loading Factor)');
ylim([-2 4]);
