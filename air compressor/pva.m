addpath('./functions/pva')

HB = 39; alpha = 90*pi/180; beta = 90*pi/180; n1 = 500; iter = 20000;
AB = 104.395; CD = 94.8311; HD = 39.0074;
OA = HB/2; AC = OA/2;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;

rAs  = zeros(iter, 1); vAs  = zeros(iter, 1); aAs  = zeros(iter, 1);
rBs  = zeros(iter, 1); vBs  = zeros(iter, 1); aBs  = zeros(iter, 1);
rCs  = zeros(iter, 1); vCs  = zeros(iter, 1); aCs  = zeros(iter, 1);
rDs  = zeros(iter, 1); vDs  = zeros(iter, 1); aDs  = zeros(iter, 1);
vc2s = zeros(iter, 1); ac2s = zeros(iter, 1); alp2s = zeros(iter, 1);
vc4s = zeros(iter, 1); ac4s = zeros(iter, 1); alp4s = zeros(iter, 1);

phis = linspace(0,360,iter);
for i = 1:iter
    rA = OA*exp(1i*phis(i)*pi/180);
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
    rD = OD*exp(1i*xOD); rDs(i) = rD; rCs(i) = rC;
    
    vA = vel(rA, 0, omg1);
    [v3_T, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
    vB = vel(rB, v3_T, 0);
    [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
    vC = vel(rC, vC_T, omgC);
    [v5_T, omg4] = vCPA2(-vC, rD, rD-rC, 0, 0);
    vD = vel(rD, v5_T, 0);
    
    aA = acc(rA, 0, omg1, 0, 0);
    [a3_T, alp2] = aCPA2(-aA, rB, rB-rA, v3_T, 0, 0, omg2, 0, 0);
    aB = acc(rB, v3_T, 0, a3_T, 0);
    [aC_T, alpC] = aCPA1(-aA, rC, rC-rA, vC_T, omgC, 0, omg2, 0, alp2);
    aC = acc(rC, vC_T, omgC, aC_T, alpC);
    [a5_T, alp4] = aCPA2(-aC, rD, rD-rC, v5_T, 0, 0, omg4, 0, 0);
    aD = acc(rD, v5_T, 0, a5_T, 0);
    
%     ----- save positions of A, B, C, D -----
    rAs(i) = rA; rBs(i) = rB; rCs(i) = rC; rDs(i) = rD;
%     sys(i, 1, :) = [link1, link2_AB, link2_AC, link2_BC, link4];

%     ----- save velocities of A, B, C, D -----
    vAs(i) = vA; vBs(i) = vB; vCs(i) = vC; vDs(i) = vD;

%     ----- save accelerations of A, B, C, D -----
    aAs(i) = aA; aBs(i) = aB; aCs(i) = aC; aDs(i) = aD;
    alp2s(i) = alp2; alp4s(i) = alp4;
end

rs = [rAs, rBs, rCs, rDs];
vs = [vAs, vBs, vCs, vDs];
as = [aAs, aBs, aCs, aDs];
writematrix([real(rs), imag(rs)], './outputs/lib_pos');
writematrix([real(vs), imag(vs)], './outputs/lib_vel');
writematrix([real(as), imag(as)], './outputs/lib_acc');
writematrix([alp2s, alp4s], './outputs/lib_ang_acc');

% ----- plot velocity + acceleration -----
% figure, hold on, xlabel('\phi (deg)'), ylabel('velocity (mm/s)')
% plot(phis, vBs, 'b')
% plot(phis, vDs, 'r'), legend('|v_B|', '|v_D|')
% figure, hold on, xlabel('\phi (deg)'), ylabel('acceleration (mm/s^2)')
% plot(phis, aBs, 'b')
% plot(phis, aDs, 'r'), legend('|a_B|', '|a_D|')

% % -----  find minimum CD ------
% rA_T = OA*exp(1i*xOB);
% [~,~,~,OB] = CPA2(xOB, AB, rA_T);
% rB_T = OB*exp(1i*xOB);
% rD_T = -real(rB_T) + 1i*imag(rB_T);
% CD = min(abs(rD_T-rCs)); % the script gives CD = 94.8311
% 
% % ----- find HD -----
% HD = max(abs(rDs)) - min(abs(rDs)); % the script gives HD = 39.0074