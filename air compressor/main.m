params = readtable('parameters.txt');
disp(params)
alpha = deg2rad(params{1,1});
beta = deg2rad(params{1,2});
n1 = params{1,3}; omg1 = convangvel(n1, 'rpm', 'rad/s');
HB = params{1,4}*1e-3;      OA = HB/2;      AC = OA;
AB = params{1,5}*1e-3;      BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
delta = params{1,17};   muE = 1;    muJ = 1;
A_piston = params{1,18}*1e-6;

formulate(alpha, OA, AC, AB, BC);

iter = 5000;
J = pvaJ(linspace(0,2*pi,iter));

%%
% For FB, Ad
sP_B = readtable('./outputs/sP_B.csv');
sB = sP_B{:,1}*1e-3;
PB = sP_B{:,2};
phi1 = radk2pi(phi1wrtPB(sB));
Md = Mdwrtphi(phi1, PB);
FB = PB*A_piston;
FB_2T  = repmat(FB,2,1);
phi_2T = [phi1;phi1+4*pi];
Md_2T  = repmat(Md,2,1);

phi = linspace(0,4*pi,720)';

[phi_u, index] = unique(phi_2T - 4*pi);
Md = interp1(phi_u, Md_2T(index), phi,'spline');
[phi_u, index] = unique(phi_2T - 4*pi); 
FB = interp1(phi_u, FB_2T(index), phi, 'spline');

Ad = cumtrapz(phi(2),Md);
phi_Ad = [phi, Ad];
FB_phi = [phi, FB];

phi_Js = [linspace(0,4*pi,iter*2)', [J;J]];

% Find Ac
sP_D = readtable('./outputs/sP_D.csv');
sD = sP_D{:,1}*1e-3;
PD = sP_D{:,2};
phi2 = radk2pi(phi1wrtPD(sD));
Mc = Mcwrtphi(phi2, PD);
rad2_2T = [phi2; phi2+2*pi; phi2+4*pi; phi2+6*pi];
Mc_2T = repmat(Mc,4,1);

[phi_u, index] = unique(rad2_2T-2*pi);
Mc = interp1(phi_u, Mc_2T(index), phi, 'spline');
Ac = cumtrapz(phi(2), Mc);

%% Matching dimensions of FD(phi), Mc(phi), J(phi) with FB(phi)
% since length(phi2) != length(phi), we match FD(phi2) with FD(phi) using
% interpolation.

FD = Ad(end)/Ac(end) * PD * A_piston;
FD_2T = repmat(FD,4,1);
FD = interp1(phi_u, FD_2T(index), phi, 'spline');

Mcr = Ad(end)/Ac(end) * Mc;
Acr = cumtrapz(phi(2), Mcr);
phi_Ac = [phi, Acr];
FD_phi = [phi, FD];

%% Find Jd using wittenbauer E-J diagram
% Find E
[phi_Js_u, index] = unique(phi_Js(:,1)); 
J = interp1(phi_Js_u, phi_Js(index,2), phi, 'spline');
E = phi_Ad(:,2) - phi_Ac(:,2) + 1/2*omg1^2*phi_Js(1,2);

% find the tangent lines y = ax + b
psi_max= atan(muJ/(2*muE)*omg1^2*(1+delta/2)^2);
psi_min= atan(muJ/(2*muE)*omg1^2*(1-delta/2)^2);
pointA=E-tan(psi_max)*J;
pointB=E-tan(psi_min)*J;
a=max(pointA);
b=min(pointB);
x_td = linspace(-0.1,0.15,10);
ytmax=tan(psi_max)*x_td+a;
ytmin=tan(psi_min)*x_td+b;
Jd=muJ*(a-b/3)/(tan(psi_max)-tan(psi_min))+phi_Js(1,2);

%% Plotting diagrams

% % FOR PB_s and PD_s
% plot(sB,P_B);
% plot(sD,PofD);
% 
% % FOR FB and FD
% plot(phi_unique,FB_unique);
% plot(phi_unique,FD_unique);
% 
% % FOR MOMENT DONG, MOMENT CAN
% plot(phi,Md);
% plot(phi,Mc);
% 
% % FOR CONG DONG, CONG CAN
% plot(phi_Ad(:,1),phi_Ad(:,2),'b',phi_Ac(:,1),phi_Ac(:,2),'r');
% 
% % FOR E-J and flywheel
% plot(phi,E,'k');

grid on, hold on
plot(J,E,'k');
plot(x_td,ytmax,'r',x_td,ytmin,'r')