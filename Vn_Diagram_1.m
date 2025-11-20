% V-n diagram
close all; clear all; clc;

g = 9.81;
c = 7.08/12; % feet
WS = 3.5;                 % lb/ft^2
rho = 0.00230;    % slugs/ft^3
Cl_max = 1.3642;
Cn_max = Cl_max*1.1;
CL3 = 4.15;
Kc = 33;                  
Ude = 12;
Ude2 = 6;                 % fps
V = linspace(0,300, 100); 

% Kg and mug
mug = 2 * (WS / (rho * c * g * CL3));
Kg = 0.88 * (mug / (5.3 + mug));

% nlims
nlimpos = 3;
nlimneg = -0.4 * nlimpos;

% nults
nultpos = 1.5 * nlimpos;
nultneg = 1.5 * nlimneg;

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
plot(V, n, 'c', 'LineWidth', 1);
plot(V, n_neg, 'c','LineWidth', 1);
yline(0, '--k', 'LineWidth', 1);
yline(nlimpos, '--', 'LineWidth', 1);
yline(nlimneg, '--', 'LineWidth', 1);
yline(nultpos, '--m', 'LineWidth', 1);
yline(nultneg, '--m', 'LineWidth', 1);
plot(V, nlimgustC, '--g', 'LineWidth', 1);
plot(V, nlimgustCneg, '--g', 'LineWidth', 1);
plot(V, nlimgustD, '--b', 'LineWidth', 1);
plot(V, nlimgustDneg, '--b', 'LineWidth', 1);
plot(Vcruise_min, 0, 'ro', 'MarkerFaceColor','r');
text(Vcruise_min, 0, '  VC min', 'HorizontalAlignment','right', 'VerticalAlignment','top', 'Color','r', 'FontWeight','bold');
plot(Vstall, 0, 'ko', 'MarkerFaceColor', 'k');
text(Vstall, 0, '  Vstall', 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'Color','k', 'FontWeight','bold');
plot(Vdive_min, 0, 'bo', 'MarkerFaceColor', 'b');
text(Vdive_min, 0, '  VD min', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Color','b', 'FontWeight','bold');
plot(Va, 0, 'ko', 'MarkerFaceColor', 'k');
text(Va, 0, '  VA', 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'Color','k', 'FontWeight','bold');
title('V-n Diagram');
xlabel('V (ft/s)');
ylabel('n (Loading Factor)');
ylim([-2 4]);
text(V(end), n(end), ' Stall +', 'Color','b', 'HorizontalAlignment','right');
text(V(end), n_neg(end), ' Stall -', 'Color','b', 'HorizontalAlignment','right');
text(V(end), nlimgustC(end), ' + VC Gust', 'Color','g', 'HorizontalAlignment','right');
text(V(end), nlimgustCneg(end), ' - VC Gust', 'Color','g', 'HorizontalAlignment','right');
text(V(end), nlimgustD(end), ' + VD Gust', 'Color','b', 'HorizontalAlignment','right');
text(V(end), nlimgustDneg(end), ' - VD Gust', 'Color','b', 'HorizontalAlignment','right');
text(10, nlimpos, ' nlim positive', 'Color','k');
text(10, nlimneg, ' nlim negative', 'Color','k');
text(10, nultpos, ' nult positive', 'Color','m');
text(10, nultneg, ' nult negative', 'Color','m');


% intersection points
diff2 = nlimgustCneg - nultneg;
idx2 = find(diff2 .* circshift(diff2,-1) <= 0, 1);   % where signs change
V_int2 = V(idx2);
n_int2 = nlimgustCneg(idx2);
xline(V_int2, 'm--', 'LineWidth', 1.5);
plot(V_int2, 0, 'ro', 'MarkerFaceColor', 'r');
text(V_int2, 0, '  VC max', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Color','r', 'FontWeight','bold');


diff3 = nlimgustDneg - nultneg;
idx3 = find(diff3 .* circshift(diff3,-1) <= 0, 1);   % where signs change
V_int3 = V(idx3);
n_int3 = nlimgustDneg(idx3);
xline(V_int3, 'm--', 'LineWidth', 1.5);
plot(V_int3, 0, 'bo', 'MarkerFaceColor', 'b');
text(V_int3, 0, '  VD max', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Color','b', 'FontWeight','bold');
