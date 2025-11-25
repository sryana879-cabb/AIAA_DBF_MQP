% V-n diagram
close all; clear all; clc;
g = 32.2;
c = 0.86985; % feet
WS = 2.6;                 % lb/ft^2
rho = 0.0023; % slugs/ft^3
Cl_max = 0.86256;
Cn_max = Cl_max*1.1;
CL3 = 0.0589 * 180/pi; % rad
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
nultpos = 1.5*nlimpos;
nultneg = 1.5*nlimneg;
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

% vertical lines at vc and vd max
diff2 = nlimgustCneg - nultneg;
idx2 = find(diff2 .* circshift(diff2,-1) <= 0, 1);   % where signs change
V_int2 = V(idx2);
n_int2 = nlimgustCneg(idx2);
diff3 = nlimgustDneg - nultneg;
idx3 = find(diff3 .* circshift(diff3,-1) <= 0, 1);   % where signs change
V_int3 = V(idx3);
n_int3 = nlimgustDneg(idx3);

% Plot
figure; hold on; grid on;

% green region
Vg = V(V > Vstall & V <= V_int2);  
Vg = [Vstall Vg];                  
n_upper = min(interp1(V, n, Vg), nlimpos);
n_lower = max(interp1(V, n_neg, Vg), nlimneg);
fill([Vg fliplr(Vg)], [n_upper fliplr(n_lower)], 'g', 'FaceAlpha',0.25, 'EdgeColor','none');

% red region
Vr = V(V >= V_int3);                        
n_red_top = nultpos * ones(size(Vr));       
n_red_bottom = nultneg * ones(size(Vr));    
fill([Vr fliplr(Vr)], [n_red_top fliplr(n_red_bottom)], 'r', 'FaceAlpha',0.25, 'EdgeColor','none');

% orange region top
V_or = V(V>=Vstall);
stall_pos_or = n(V>=Vstall);
upper_or = max(stall_pos_or, nlimpos);
fill([V_or fliplr(V_or)], [upper_or fliplr(nlimpos*ones(size(upper_or)))], [1 0.5 0], 'FaceAlpha',0.30,'EdgeColor','none');

% orange region bottom
Vo_bot = V(V >= Vstall & V <= V_int3);               % same V range
n_orange_bottom2 = max(interp1(V, n_neg, Vo_bot), nultneg);  % top boundary = nultneg
n_orange_top2 = nlimneg * ones(size(Vo_bot));       % bottom boundary = nlimneg
fill([Vo_bot fliplr(Vo_bot)], [n_orange_top2 fliplr(n_orange_bottom2)], [1 0.5 0], 'FaceAlpha',0.30, 'EdgeColor','none');

% yellow region
Vy = V(V >= V_int2 & V <= V_int3);
y_upper = nlimpos * ones(size(Vy));
y_lower = nlimneg * ones(size(Vy));
fill([Vy fliplr(Vy)], [y_upper fliplr(y_lower)], 'y', 'FaceAlpha',0.25, 'EdgeColor','none');

% plot lines
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
text(Va, 0, '  VA', 'HorizontalAlignment','right', 'VerticalAlignment','top', 'Color','k', 'FontWeight','bold');
title('V-n Diagram');
xlabel('V (ft/s)');
ylabel('n (Loading Factor)');
ylim([-2 4]);
text(V(end), n(end), ' Stall +', 'Color','b', 'HorizontalAlignment','right');
text(V(end), n_neg(end), ' Stall -', 'Color','b', 'HorizontalAlignment','right');
text(V(10), nlimgustC(10), ' + VC Gust', 'Color','g', 'HorizontalAlignment','right');
text(V(10), nlimgustCneg(10), ' - VC Gust', 'Color','g', 'HorizontalAlignment','right');
text(V(10), nlimgustD(10), ' + VD Gust', 'Color','b', 'HorizontalAlignment','right');
text(V(10), nlimgustDneg(10), ' - VD Gust', 'Color','b', 'HorizontalAlignment','right');
text(10, nlimpos, ' nlim positive', 'Color','k');
text(10, nlimneg, ' nlim negative', 'Color','k');
text(10, nultpos, ' nult positive', 'Color','m');
text(10, nultneg, ' nult negative', 'Color','m');
% intersection points
plot(V_int2, 0, 'ro', 'MarkerFaceColor', 'r');
text(V_int2, 0, '  VC max', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Color','r', 'FontWeight','bold');
xline(V_int2, 'm--', 'LineWidth', 1.5);
plot(V_int3, 0, 'bo', 'MarkerFaceColor', 'b');
text(V_int3, 0, '  VD max', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Color','b', 'FontWeight','bold');
xline(V_int3, 'm--', 'LineWidth', 1.5);



