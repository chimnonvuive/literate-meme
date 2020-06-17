function [Js, rAs, rBs, rCs, rDs] = pvaJ(phis)

    params = readtable('parameters');
    alpha = params{1,1}*pi/180;     xOB = pi/2 - alpha/2;
    beta  = params{1,2}*pi/180;     xOD = pi/2 + alpha/2; 
    n1 = params{1,3}; omg1 = n1*pi/30;
    HB = params{1,4}*1e-3;      OA = HB/2;      AC = OA;
    AB = params{1,5}*1e-3;      BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
    CD = params{1,6}*1e-3;

    m1 = params{1,7};         J1 = params{1,13};
    m2 = params{1,8};         J2 = params{1,14};
    m3 = params{1,9};         J3 = params{1,15};
    m4 = params{1,10};        J4 = params{1,16};
    m5 = params{1,11};        J5 = params{1,17};
    
    iter = length(phis);
    rAs  = zeros(iter, 1); vAs  = zeros(iter, 1); aAs  = zeros(iter, 1);
    rBs  = zeros(iter, 1); vBs  = zeros(iter, 1); aBs  = zeros(iter, 1);
    rCs  = zeros(iter, 1); vCs  = zeros(iter, 1); aCs  = zeros(iter, 1);
    rDs  = zeros(iter, 1); vDs  = zeros(iter, 1); aDs  = zeros(iter, 1);
    alp2s = zeros(iter, 1); Js = zeros(iter,1);
    alp4s = zeros(iter, 1);
    
    for i = 1:iter
        rA = OA*exp(1i*phis(i));
        [~,~,~,OB] = CPA2(xOB, AB, rA);
        rB = OB*exp(1i*xOB);
        [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
        rC = rA + AC*exp(1i*xAC);
        [~,~,~,OD] = CPA2(xOD, CD, rC);
        rD = OD*exp(1i*xOD);
         rDs(i) = rD;
         rCs(i) = rC;

        vA = vel(rA, 0, omg1);
        [v3_T, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
        vB = vel(rB, v3_T, 0);
        [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
        vC = vel(rC, vC_T, omgC);
        [v5_T, omg4] = vCPA2(-vC, rD, rD-rC, 0, 0);
        vD = vel(rD, v5_T, 0);
        
        centerCA = (rC+rA)/2;
        centerAB = rA*0/OA; % Use Solidworks to find COM
        rS2 = centerCA+(centerAB-centerCA)*AC/AB;
        rS4 = (rC+rD)/2;   

        vS2 = vel(rS2, 0, omg2);
        vS4 = vel(rS4, 0, omg4);
        
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
        
        Js(i) = J1 + J2*(omg2/omg1)^2 + m2*(abs(vS2)/omg1)^2 + ...
            m3*(abs(vB)/omg1)^2 + J4*(omg4/omg1)^2 + ...
            m4*(abs(vS4)/omg1)^2 + m5*(abs(vD)/omg1)^2;
    end
    rs = [rAs, rBs, rCs, rDs];
    vs = [vAs, vBs, vCs, vDs];
    as = [aAs, aBs, aCs, aDs];
    writematrix([real(rs), imag(rs)], './outputs/lib_pos.txt');
    writematrix([real(vs), imag(vs)], './outputs/lib_vel.txt');
    writematrix([real(as), imag(as)], './outputs/lib_acc.txt');
    writematrix([alp2s, alp4s], './outputs/lib_ang_acc.txt');
    writematrix(Js, './outputs/lib_J.txt');
end

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