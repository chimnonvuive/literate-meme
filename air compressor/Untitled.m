load('F_D_and_phi.mat', 'F_D_and_phi');
P_Ds = F_D_and_phi(:,2); phis = F_D_and_phi(:,1);
P_Bs = readtable('FB_new.txt'); P_Bs = P_Bs{:,1};
A_piston = 100; % cm^2
F_Ds = P_Ds*A_piston/1e4; % (N)

% ----- Gravity on 5 links -----
m1 = 1;    m2 = 1;    m3 = 1;    m4 = 1;    m5 = 1;    g = -9.81j;
Q1 = m1*g; Q2 = m2*g; Q3 = m3*g; Q4 = m4*g; Q5 = m5*g;

% ----- Designed parameters -----
HB = 39; HD = 39.0074; OA = HB/2; AB = 104.395; AC = OA/2; CD = 94.8311;
alpha = 90; beta = 90; n1 = 500;
BC = sqrt(AC^2+AB^2-2*AC*AB*cosd(beta));
xOB = deg2rad(90 - alpha/2); xOD = deg2rad(90 + alpha/2);
omg1 = convangvel(n1,'rpm','rad/s');

iter = length(F_Ds);
%% SECTION TITLE
% DESCRIPTIVE TEXT



% new_Ad= zeros(length(Ac_energy),1);
% [x, index] = unique(Ad_energy(:,1));
% new_Ad=interp1(x,Ad_energy(index,2),Ac_energy(:,1));
% F_B_wrt_phiD = interp1(phiBs,F_Bs,phiDs);
