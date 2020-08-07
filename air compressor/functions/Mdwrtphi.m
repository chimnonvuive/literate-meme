function Mds = Mdwrtphi(phis, PB)
    
    params = readtable('parameters');
    CD = readtable('./outputs/CD.txt');
    load('./outputs/functions.mat');
    
    alpha = deg2rad(params{1,1});     xOB = pi/2 - alpha/2;
    omg1 = convangvel(params{1,3}, 'rpm', 'rad/s');
    HB = params{1,4}*1e-3;      OA = HB/2;      AC = OA;
    AB = params{1,5}*1e-3;      CD = CD{1,1};
    
    g  = params{1,11};
    m1 = params{1,6};   Q1 = -m1*g*[0,1,0];
    m2 = params{1,7};   Q2 = -m2*g*[0,1,0];
    m3 = params{1,8};   Q3 = -m3*g*[0,1,0];
    m4 = params{1,9};   Q4 = -m4*g*[0,1,0];
    m5 = params{1,10};  Q5 = -m5*g*[0,1,0];
    
    A_piston = params{1,18}*1e-6;
    mAC = AC/(AB+AC)*m2;
    mAB = AB/(AB+AC)*m2;
    
    iter = length(phis);
    Mds = zeros(iter, 1);
    
    for i=1:iter
        
        F_B = PB(i)*A_piston*[cos(xOB+pi), sin(xOB+pi), 0];
        vA = fvA(phis(i));
        vB = fvB(phis(i));
        vC = fvC(phis(i));
        vD = fvD(phis(i), CD);
        
        vAC = (vA + vC)/2;
        vAB = (vA + vB)/2;
        
        v1 = vA/2;
        v2 = (mAC*vAC + mAB*vAB)/m2;
        v4 = (vC + vD)/2;

        Mds(i) = (dot(F_B, vB) + dot(Q1, v1) + dot(Q2, v2) + dot(Q3, vB) +...
            dot(Q4, v4) + dot(Q5, vD))/omg1;
    end
end