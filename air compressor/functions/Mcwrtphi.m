function [Mcs, F_Ds] = Mcwrtphi(phi1, PD)
    
    params = readtable('parameters');
    alpha = deg2rad(params{1,1});     xOB = pi/2 - alpha/2;
    beta  = deg2rad(params{1,2});     xOD = pi/2 + alpha/2; 
    omg1 = convangvel(params{1,3}, 'rpm', 'rad/s');
    HB = params{1,4}*1e-3;      OA = HB/2;      AC = OA;
    AB = params{1,5}*1e-3;      BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
    CD = params{1,6}*1e-3;
    d_piston = params{1,19}*1e-3;   A_piston = d_piston^2*pi/4;

    iter = length(phi1);
    Mcs = zeros(iter, 1);
    F_Ds = PD*A_piston*exp(1j*(xOD+pi));
    
    for i=1:iter
        rA = OA*exp(1j*phi1(i));
        [~,~,~,OB] = CPA2(xOB, AB, rA);
        rB = OB*exp(1j*xOB);
        [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
        rC = rA + AC*exp(1j*xAC);
        [~,~,~,OD] = CPA2(xOD, CD, rC);
        rD = OD*exp(1j*xOD);

        vA = vel(rA, 0, omg1);
        [~, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
        [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
        vC = vel(rC, vC_T, omgC);
        [v5_T, ~] = vCPA2(-vC, rD, rD-rC, 0, 0);
        vD = vel(rD, v5_T, 0);
        
        Mcs(i)=cdot(F_Ds(i),vD)/omg1;
    end
end