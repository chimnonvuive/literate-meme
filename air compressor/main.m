params = readtable('parameters.txt');
disp(params)
alpha = deg2rad(params{1,1});
beta = deg2rad(params{1,2});
n1 = params{1,3}; omg1 = convangvel(n1, 'rpm', 'rad/s');
HB = params{1,4}*1e-3;
OA = HB/2;
AC = params{1,5}*1e-3;
AB = 5.372*OA;
BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
delta = params{1,17};   muE = 1;    muJ = 1;
A_piston = params{1,18}*1e-6;

% the center of gravity of link 1 is located at point O
% the center of gravity of link 2 is located at the midpoint of AB

formulate(omg1, alpha, OA, AC, AB, BC);

iter = 500;
phi = linspace(0, 2*pi, iter)';
J = pvaJ(phi);

%%
% Find Md(phi), FB(phi) and Ad(phi) at B
sP_B = readtable('./outputs/sP_B.csv');
sB = sP_B{:,1}*1e-3;
PB = sP_B{:,2};
% angles of link 1 in 1 periods corresponding to sB
phi1 = radk2pi(phi1wrtPB(sB, OA, AB, alpha));
% pressure at B in 2 periods
FB  = repmat(PB*A_piston,2,1);
% dynamic moment in 2 periods
Md  = repmat(Mdwrtphi(phi1, PB), 2, 1);

% angles of link 1 with equal angular displacement in 2 periods
phi = [phi; phi+2*pi];
% moment of inertia in 2 periods
J = [J; J];
% the corresponding Md(phi), FB(phi) and Ad(phi)
[phi_u, index] = unique([phi1-4*pi; phi1]);
Md = interp1(phi_u, Md(index), phi, 'spline');
FB = interp1(phi_u, FB(index), phi, 'spline');
Ad = cumtrapz(phi(2), Md); % Ad = Md*phi

%% Find Mc(phi), FD(phi) and Ac(phi) at D
sP_D = readtable('./outputs/sP_D.csv');
sD = sP_D{:,1}*1e-3;
PD = sP_D{:,2};
% angles of link 1 in 1 periods corresponding to sD
phi1 = radk2pi(phi1wrtPD(sD, OA, AB, alpha));
% pressure at D in 4 half periods
Mc = repmat(Mcwrtphi(phi1, PD), 4, 1);

% generate data for 4 half periods
[phi_u, index] = unique([phi1-2*pi; phi1; phi1+2*pi; phi1+4*pi]);
Mc = interp1(phi_u, Mc(index), phi, 'spline');
Ac = cumtrapz(phi(2), Mc); % Ac = Mc*phi
% pressure at D in 4 half periods
FD = repmat(Ad(end)/Ac(end) * -PD * A_piston, 4, 1);

% the corresponding Mc(phi), FD(phi) and Ac(phi)
Mc = Mc * Ad(end)/Ac(end);
FD = interp1(phi_u, FD(index), phi, 'spline');
Ac = cumtrapz(phi(2), Mc);

% Find E
[phi_u, index] = unique(phi); 
J = interp1(phi_u, J(index), phi, 'spline');
E = Ad - Ac + 1/2*omg1^2*J(1);

% find the tangent lines y = ax + b
psi_max = atan( muJ/(2*muE) * omg1^2 * (1+delta/2)^2 );
psi_min = atan( muJ/(2*muE) * omg1^2 * (1-delta/2)^2 );
a = max(E - tan(psi_max)*J);
b = min(E - tan(psi_min)*J);
x = linspace(-10,10,100);
ymax = tan(psi_max)*x + a;
ymin = tan(psi_min)*x + b;
Jd = muJ*(a-b)/(tan(psi_max)-tan(psi_min)) + J(2); disp(Jd)
% clear a b A B psi_max psi_min phi1 phi_u index

%% Plotting diagrams

% figure, hold on, grid on
% title('Displacement - pressure relation at B')
% xlabel('displacement s_B (mm)'), ylabel('pressure P_B (bar)')
% plot(sB*1e3, PB*1e-5, 'k-', 'linewidth', 1.5), legend('pressure input')
% hold off
% 
% figure, hold on, grid on 
% title('Displacement - pressure relation at D')
% xlabel('displacement s_D (mm)'), ylabel('pressure P_D (bar)')
% plot(sD*1e3, PD*1e-5, 'k-', 'linewidth', 1.5), legend('resistant pressure')
% hold off

figure, hold on, grid on
title('Force applied on 2 pistons in 1 cycle')
xlabel('angle \phi_1 (deg)'), ylabel('force (kN)')
xlim([0, 720]), xticks(0:90:720)
plot(rad2deg(phi), FB*1e-3, 'b-', 'linewidth', 1.5)
plot(rad2deg(phi), FD*1e-3, 'r-', 'linewidth', 1.5)
legend('piston B', 'piston D')
hold off

figure, hold on, grid on
title('Moment applied on 2 pistons in 1 cycle')
xlabel('angle \phi_1 (deg)'), ylabel('moment (N\cdotmm)')
xlim([0, 720]), xticks(0:90:720)
plot(rad2deg(phi), Md*1e3, 'b-', 'linewidth', 1.5)
plot(rad2deg(phi), Mc*1e3, 'r-', 'linewidth', 1.5)
legend('dynamic moment on piston B', 'resistant moment on piston D')
hold off

figure, hold on, grid on
title('Work applied on 2 pistons in 1 cycle')
xlabel('angle \phi_1 (deg)'), ylabel('work (J)')
xlim([0, 720]), xticks(0:90:720)
plot(rad2deg(phi), Ad, 'b-', 'linewidth', 1.5)
plot(rad2deg(phi), Ac, 'r-', 'linewidth', 1.5)
legend('dynamic work on piston B', 'resistant work on piston D')
hold off 

figure, hold on, grid on
title('Energy output at D in 1 cycle')
xlabel('angle \phi_1 (deg)'), ylabel('Energy (J)')
xlim([0, 720]), xticks(0:90:720)
plot(rad2deg(phi), E, 'k-', 'linewidth', 1.5)
legend('energy output')
hold off

figure, hold on, grid on
title('Wittenbauer graph of the system in 1 cycle')
xlabel('J (kg\cdotm^2)'), ylabel('E (J)')
xlim([min(J)*0.9, max(J)*1.1]), ylim([min(E)*0.9, max(E)*1.05])
plot(J, E, 'k-', 'linewidth', 1.5)
plot(x, ymax, 'r', x, ymin, 'r')
legend('Wittenbauer curve', 'tangential boundaries')
hold off
%%
vD = readtable('./outputs/lib_acc.txt');
figure, hold on, grid on
title('Acceleration along x-axis of piston B')
xlim([0, 360]), xticks(0:90:360)
xlabel('angle \phi_1 (deg)'), ylabel('a_{Bx} (mm/s^2)')
plot(linspace(0,360,length(vD{:,2})), vD{:,2}*1e3, 'k-', 'linewidth', 1.5)
legend('piston B')
hold off