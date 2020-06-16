%% Generated external forces and corresponding pva
data = readtable('PD_phiD.txt'); phis = data{:,1}; P_Ds = data{:,2};
P_Bs = readtable('PB_wrt_phiD.txt'); P_Bs = P_Bs{:,1}*10;

rs = readtable('system_positions.txt');
rAs = rs{:,1} + rs{:,5}*1j;
rBs = rs{:,2} + rs{:,6}*1j;
rCs = rs{:,3} + rs{:,7}*1j;
rDs = rs{:,4} + rs{:,8}*1j;

vs = readtable('system_velocities.txt');
vAs = vs{:,1} + vs{:,5}*1j;
vBs = vs{:,2} + vs{:,6}*1j;
vCs = vs{:,3} + vs{:,7}*1j;
vDs = vs{:,4} + vs{:,8}*1j;

as = readtable('system_accelerations.txt');
aAs = as{:,1} + as{:,5}*1j;
aBs = as{:,2} + as{:,6}*1j;
aCs = as{:,3} + as{:,7}*1j;
aDs = as{:,4} + as{:,8}*1j;

alps = readtable('system_angular_accs.txt');
alp2s = alps{:,1}; alp4s = alps{:,2};
%% ----- Gravity on 5 links -----
m1 = 1;    m2 = 1;    m3 = 1;    m4 = 1;    m5 = 1;    g = -9.81j;
Q1 = m1*g; Q2 = m2*g; Q3 = m3*g; Q4 = m4*g; Q5 = m5*g;

%% ----- Design parameters -----
HB = 39; HD = 39.0074; OA = HB/2; AB = 104.395; AC = OA/2; CD = 94.8311;
alpha = 90; beta = 90; n1 = 500; A_piston_cm = 100; I2 = 0.12; I4 = 0.12;
BC = sqrt(AC^2+AB^2-2*AC*AB*cosd(beta));
xOB = 90 - alpha/2; xOD = 90 + alpha/2;
omg1 = convangvel(n1,'rpm','rad/s');

F_Bs = P_Bs*A_piston_cm/1e4*exp(1j*-deg2rad(xOB)); % (N)
F_Ds = P_Ds*A_piston_cm/1e4*exp(1j*-deg2rad(xOD)); % (N)
Fs_eq_Ms = zeros(length(phis), 8);

for i=1:length(phis)
    a4 = (aDs(i)+aCs(i))/2; a2 = (aAs(i)+aBs(i)+aCs(i))/3;
    rAG1 = rAs(i) - (rAs(i)+0)/2;
    rOG1 = 0 - (rAs(i)+0)/2;
    rDG4 = rDs(i) - (rDs(i)+rCs(i))/2;
    rCG4 = rCs(i) - (rDs(i)+rCs(i))/2;
    rAG2 = rAs(i) - (rAs(i)+rBs(i)+rCs(i))/3;
    rBG2 = rBs(i) - (rAs(i)+rBs(i)+rCs(i))/3;
    rCG2 = rCs(i) - (rAs(i)+rBs(i)+rCs(i))/3;
%        re(F05)  im(F05)     re(F45)      im(F45)      re(F24)     im(F24)
   link45 = [ 1,  0,          1,           0,           0,          0;
              0,  1,          0,           1,           0,          0;
             -1,  1,          0,           0,           0,          0;
              0,  0,         -1,           0,           1,          0;
              0,  0,          0,          -1,           0,          1;
              0,  0, imag(rDG4), -real(rDG4), -imag(rCG4), real(rCG4)];
   link45b = [real(m5*aDs(i) - F_Ds(i) - Q5);
              imag(m5*aDs(i) - F_Ds(i) - Q5);
              0;
              real(m4*a4 - Q4);
              imag(m4*a4 - Q4);
              I4*alp4s(i)];
   res = linsolve(link45, link45b);
   F05 = res(1) + 1j*res(2);
   F45 = res(3) + 1j*res(4);
   F24 = res(5) + 1j*res(6);
   
%        re(F12)  im(F12)     re(F23)      im(F23)      re(F03)     im(F03)
   link23 = [ 0,  0,          1,           0,           1,          0;
              0,  0,          0,           1,           0,          1;
              0,  0,          0,           0,           1,          1;
              -imag(rAG2), real(rAG2), imag(rBG2), -real(rBG2), 0,  0;
              1,  0,         -1,           0,           0,          0;
              0,  1,          0,          -1,           0,          0];
   link23b = [real(m3*aBs(i) - F_Bs(i) - Q3);
              imag(m3*aBs(i) - F_Bs(i) - Q3);
              0;
              I2*alp2s(i) - (-imag(rCG2)*real(-F24)+real(rCG2)*imag(-F24));
              real(m2*a2 + F24 - Q2);
              imag(m2*a2 + F24 - Q2)];
   res = linsolve(link23, link23b);
   F12 = res(1) + 1j*res(2);
   F23 = res(3) + 1j*res(4);
   F03 = res(5) + 1j*res(6);
   F01 = m1*aAs(i)/2 - -F12 - Q1;
   M_eq = -(-imag(rAG1)*real(-F12)+real(rAG1)*imag(-F12) + -imag(rOG1)*real(F01)+real(rOG1)*imag(F01));
   
   Fs_eq_Ms(i,:) = [F01, F03, F05, F12, F23, F24, F45, M_eq];
end