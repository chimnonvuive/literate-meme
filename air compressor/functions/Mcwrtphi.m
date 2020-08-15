function Mcs = Mcwrtphi(phis, PD)
    
    params = readtable('parameters');
    CD = readtable('./outputs/length_cd.txt'); CD = CD{1,1};
    load('./outputs/functions.mat');
    
    alpha = deg2rad(params{1,1});     xOD = pi/2 + alpha/2;
    omg1 = convangvel(params{1,3}, 'rpm', 'rad/s');
    A_piston = params{1,18}*1e-6;

    iter = length(phis);
    Mcs = zeros(iter, 1);
    
    for i=1:iter
        vD = fvD(phis(i), CD);
        F_D = PD(i)*A_piston*[cos(xOD+pi), sin(xOD+pi), 0];
        Mcs(i) = dot(F_D, vD)/omg1;
    end
end