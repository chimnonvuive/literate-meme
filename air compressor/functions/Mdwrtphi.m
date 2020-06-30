function [Mds, F_Bs] = Mdwrtphi(phi1, PB)
    
    params = readtable('parameters');
    alpha = deg2rad(params{1,1});     xOB = pi/2 - alpha/2;
    beta  = deg2rad(params{1,2});     xOD = pi/2 + alpha/2; 
    omg1 = convangvel(params{1,3}, 'rpm', 'rad/s');
    HB = params{1,4}*1e-3;      OA = HB/2;      AC = OA;
    AB = params{1,5}*1e-3;      BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
    CD = params{1,6}*1e-3;
    g  = params{1,12};
    m2 = params{1,8};   Q2 = -1j*m2*g;
    m3 = params{1,9};   Q3 = -1j*m3*g;
    m4 = params{1,10};  Q4 = -1j*m4*g;
    m5 = params{1,11};  Q5 = -1j*m5*g;
    d_piston = params{1,19}*1e-3;   A_piston = d_piston^2*pi/4;

    iter = length(phi1);
    Mds = zeros(iter, 1);
    F_Bs = PB*A_piston*exp(1j*(xOB+pi));
    
    for i=1:iter
        rA = OA*exp(1j*phi1(i));
        [~,~,~,OB] = CPA2(xOB, AB, rA);
        rB = OB*exp(1j*xOB);
        [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
        rC = rA + AC*exp(1j*xAC);
        [~,~,~,OD] = CPA2(xOD, CD, rC);
        rD = OD*exp(1j*xOD);

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
        Mds(i) = (cdot(F_Bs(i),vB) + cdot(Q3,vB) + cdot(Q2,vS2) + ...
            cdot(Q4,vS4) + cdot(Q5,vD))/omg1;
    end
end