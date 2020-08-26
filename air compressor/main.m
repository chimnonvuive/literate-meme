global phi
params = readtable('parameters.txt');
disp(params)
alpha = deg2rad(params{1,1});
beta = deg2rad(params{1,2});
n1 = params{1,3}; omg1 = n1*pi/30; T1 = 2*pi/omg1;
HB = params{1,4}*1e-3;
OA = HB/2;
AC = params{1,5}*1e-3;
AB = 5.372*OA;
BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
delta = params{1,17};   muE = 1;    muJ = 1;
A_piston = params{1,18}*1e-6;

% the center of gravity of link 1 is located at point O
% the center of gravity of link 2 is located at the midpoint of AB

formulate(omg1, OA, AC, AB, BC);

iter = 500;
phi = linspace(0, 2*pi, iter)';
J = pvaJ(phi);

%% 1 cycle = 2 periods
phi = linspace(0, 4*pi, iter*2)'; J = [J; J];

sP_B = readtable('./outputs/sP_B.csv');
sB = sP_B{:,1}*1e-3; PB = sP_B{:,2};
% angles of link 1 in 1 periods corresponding to sB
phi1 = radk2pi(phi1wrtPB(sB, OA, AB, alpha));
phi2T = [phi1-4*pi; phi1];
FB = wrtphi(phi2T, repmat(PB*A_piston,2,1));
Md = wrtphi(phi2T, repmat(Mdwrtphi(phi1, PB), 2, 1));
Ad = cumtrapz(phi(2), Md); % Ad = Md*phi

sP_D = readtable('./outputs/sP_D.csv');
sD = sP_D{:,1}*1e-3; PD = sP_D{:,2};
% angles of link 1 in 1 periods corresponding to sD
phi1 = radk2pi(phi1wrtPD(sD, OA, AB, alpha));
phi2T = [phi1-2*pi; phi1; phi1+2*pi; phi1+4*pi];
Mc = wrtphi(phi2T, repmat(Mcwrtphi(phi1, PD), 4, 1));
Ac = cumtrapz(phi(2), Mc); % Ac = Mc*phi
correction = Ad(end)/Ac(end);
FD = wrtphi(phi2T, repmat(correction * PD * A_piston, 4, 1));
Mc = Mc * correction;
Ac = cumtrapz(phi(2), Mc);

E = Ad - Ac + 1/2*omg1^2*J(1);

% find the tangent lines y = ax + b
psi_max = muJ/(2*muE) * omg1^2 * (1+delta/2)^2;
psi_min = muJ/(2*muE) * omg1^2 * (1-delta/2)^2;
a_max = max(E - psi_max*J);
a_min = min(E - psi_min*J);
x = linspace(-10,10,100);
ymax = psi_max*x + a_max;
ymin = psi_min*x + a_min;
Jd = muJ*(a_max-a_min)/(psi_max-psi_min) + J(2);
fprintf('flywheel moment of inertia = %f\n', Jd)
dlmwrite('./outputs/t_FB.dat',[linspace(0,2*T1,iter*2)', FB]);
dlmwrite('./outputs/t_FD.dat',[linspace(0,2*T1,iter*2)', FD]);

%% Plotting diagrams

% figure, hold on, grid on
% title('Displacement - pressure relation at B')
% xlabel('displacement s_B (mm)'), ylabel('pressure P_B (bar)')
% plot(sB*1e3, PB*1e-5, 'k-', 'linewidth', 1.5), legend('pressure input')
% hold off, print('-dpng','-r600','./outputs/sB_PB')
% 
% figure, hold on, grid on 
% title('Displacement - pressure relation at D')
% xlabel('displacement s_D (mm)'), ylabel('pressure P_D (bar)')
% plot(sD*1e3, PD*1e-5, 'k-', 'linewidth', 1.5), legend('resistant pressure')
% hold off, print('-dpng','-r600','./outputs/sD_PD')
% 
% figure, hold on, grid on
% title('Force applied on  piston B in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('force (kN)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), FB*1e-3, 'k-', 'linewidth', 1.5)
% legend('piston B')
% hold off, print('-dpng','-r600','./outputs/FB')
% 
figure, hold on, grid on
title('Force applied on piston D in 1 cycle')
xlabel('angle \phi_1 (deg)'), ylabel('force (kN)')
xlim([0, 720]), xticks(0:90:720)
plot(rad2deg(phi), FD*1e-3, 'k-', 'linewidth', 1.5)
legend('piston D')
hold off, print('-dpng','-r600','./outputs/FD')
% 
% figure, hold on, grid on
% title('Moment applied on link 1 in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('moment (N\cdotm)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), Md, 'k-', 'linewidth', 1.5)
% legend('dynamic moment')
% hold off, print('-dpng','-r600','./outputs/MB')
% 
% figure, hold on, grid on
% title('Moment applied on link 1 in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('moment (N\cdotm)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), Mc, 'k-', 'linewidth', 1.5)
% legend('resistant moment')
% hold off, print('-dpng','-r600','./outputs/MD')
% 
% figure, hold on, grid on
% title('Work applied on piston B in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('work (J)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), Ad, 'k-', 'linewidth', 1.5)
% legend('dynamic work on piston B')
% hold off, print('-dpng','-r600','./outputs/AB')
% 
% figure, hold on, grid on
% title('Work applied on piston D in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('work (J)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), Ac, 'k-', 'linewidth', 1.5)
% legend('resistant work on piston D')
% hold off, print('-dpng','-r600','./outputs/AD')
% 
% figure, hold on, grid on
% title('Work applied on 2 pistons in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('work (J)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), Ac, 'b-', 'linewidth', 1.5)
% plot(rad2deg(phi), Ad, 'r-', 'linewidth', 1.5)
% legend('resistant work on piston D', 'dynamic work on piston B')
% hold off, print('-dpng','-r600','./outputs/ABD')
% 
% figure, hold on, grid on
% title('Energy output in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('Energy (J)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), E, 'k-', 'linewidth', 1.5)
% legend('energy output')
% hold off, print('-dpng','-r600','./outputs/E')
% 
% figure, hold on, grid on
% title('Equivalent moment of inertia in 1 cycle')
% xlabel('angle \phi_1 (deg)'), ylabel('Moment of inertia (kg\cdot m^2)')
% xlim([0, 720]), xticks(0:90:720)
% plot(rad2deg(phi), J, 'k-', 'linewidth', 1.5)
% legend('energy output')
% hold off, print('-dpng','-r600','./outputs/J')
% 
% figure, hold on, grid on
% title('Energy output - equivalent moment of inertia relation')
% xlabel('J (kg\cdotm^2)'), ylabel('E (J)')
% xlim([min(J)*0.95, max(J)*1.05]), ylim([min(E), max(E)*1.05])
% plot(J, E, 'k-', 'linewidth', 1.5)
% legend('Wittenbauer curve')
% hold off, print('-dpng','-r600','./outputs/JE')
% 
% figure, hold on, grid on
% title('Energy output - equivalent moment of inertia relation')
% title('Wittenbauer graph of the system in 1 cycle')
% xlabel('J (kg\cdotm^2)'), ylabel('E (J)')
% xlim([0, max(J)*1.1]), ylim([-20, max(E)*1.1])
% plot(J, E, 'k-', 'linewidth', 1.5)
% plot(x, ymax, 'r', x, ymin, 'r')
% legend('Wittenbauer curve', 'tangential boundaries')
% hold off, print('-dpng','-r600','./outputs/JEline')
% 
% vB = readtable('./outputs/lib_vel.txt');
% figure, hold on, grid on
% title('Velocity along x-axis of piston B')
% xlim([0, 360]), xticks(0:90:360)
% xlabel('angle \phi_1 (deg)'), ylabel('v_{Bx} (mm/s)')
% plot(linspace(0,360,length(vB{:,2})), vB{:,2}*1e3, 'k-', 'linewidth', 1.5)
% legend('piston B')
% hold off, print('-dpng','-r600','./outputs/velB')
% 
% aB = readtable('./outputs/lib_acc.txt');
% figure, hold on, grid on
% title('Acceleration along x-axis of piston B')
% xlim([0, 360]), xticks(0:90:360)
% xlabel('angle \phi_1 (deg)'), ylabel('a_{Bx} (mm/s^2)')
% plot(linspace(0,360,length(aB{:,2})), aB{:,2}*1e3, 'k-', 'linewidth', 1.5)
% legend('piston B')
% hold off, print('-dpng','-r600','./outputs/accB')

clear a_max a_min iter psi_max psi_min phi1 phi2T sP_B sP_D vB aB

function out = wrtphi(phi_data, y_data)
    global phi
    [phi_u, index] = unique(phi_data);
    out = interp1(phi_u, y_data(index), phi, 'pchip');
end