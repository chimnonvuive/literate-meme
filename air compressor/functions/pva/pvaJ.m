function [Js, rAs, rBs, rCs, rDs] = pvaJ(phis)

    params = readtable('parameters');
    load('./outputs/functions.mat');
    
    n1 = params{1,3}; omg1 = convangvel(n1, 'rpm', 'rad/s');
    omg1 = [0, 0, omg1];
    HB = params{1,4}*1e-3;      OA = HB/2;

    m1 = params{1,6};         J1 = params{1,12};
    m2 = params{1,7};         J2 = params{1,13};
    m3 = params{1,8};         J3 = params{1,14};
    m4 = params{1,9};         J4 = params{1,15};
    m5 = params{1,10};        J5 = params{1,16};
    
    iter = length(phis);
    rAs  = zeros(iter, 3); vAs  = zeros(iter, 3); aAs  = zeros(iter, 3);
    rBs  = zeros(iter, 3); vBs  = zeros(iter, 3); aBs  = zeros(iter, 3);
    rCs  = zeros(iter, 3); vCs  = zeros(iter, 3); aCs  = zeros(iter, 3);
    rDs  = zeros(iter, 3); vDs  = zeros(iter, 3); aDs  = zeros(iter, 3);
    
    omg2s = zeros(iter, 3);
    omg4s = zeros(iter, 3);
    alp2s = zeros(iter, 3);
    alp4s = zeros(iter, 3);
    Js = zeros(iter,1);
    
    for i = 1:iter
        rAs(i, :) = fpA(phis(i));
        rBs(i, :) = fpB(phis(i));
        rCs(i, :) = fpC(phis(i));
    end
    
    % Find CD (at time t0 when OA is perpendicular to Ox)
    rA_perp = [0, OA, 0];
    dist = euclidean(rA_perp, rAs);
    index = find(dist == min(dist));
    rB_perp = rBs(index, :);
    rD_perp = [-rB_perp(1), rB_perp(2), 0];
    CD = min(euclidean(rD_perp, rCs));
    
    for i = 1:iter
        rDs(i, :) = fpD(phis(i), CD);
        
    %     ----- save velocities of A, B, C, D -----
        
        vAs(i, :) = fvA(phis(i));
        vBs(i, :) = fvB(phis(i));
        vCs(i, :) = fvC(phis(i));
        vDs(i, :) = fvD(phis(i), CD);
        
    %     ----- save accelerations of A, B, C, D -----
     
        aAs(i, :) = faA(phis(i));
        aBs(i, :) = faB(phis(i));
        aCs(i, :) = faC(phis(i));
        aDs(i, :) = faD(phis(i), CD);
        
        v2 = (vAs(i, :) + vBs(i, :))/2;
        v3 = vBs(i, :);
        v4 = (rCs(i, :) + rDs(i, :))/2;
        v5 = vDs(i, :);
                
        omg2s(i, :) = fomg2(phis(i)); omg2 = omg2s(i, :);
        omg4s(i, :) = fomg4(phis(i), CD);
        alp2s(i, :) = falp2(phis(i)); omg4 = omg4s(i, :);
        alp4s(i, :) = falp4(phis(i), CD);

         Js(i) = (J1*dot(omg1, omg1) +                 ...
                  J2*dot(omg2, omg2) + m2*dot(v2, v2) +...
                                       m3*dot(v3, v3) +...
                  J4*dot(omg4, omg4) + m4*dot(v4, v4) +...
                                     + m5*dot(v5, v5) )/dot(omg1, omg1);
    end
    
    rs = [rAs(:,1), rBs(:,1), rCs(:,1), rDs(:,1),...
        rAs(:,2), rBs(:,2), rCs(:,2), rDs(:,2)];
    vs = [vAs(:,1), vBs(:,1), vCs(:,1), vDs(:,1),...
        vAs(:,2), vBs(:,2), vCs(:,2), vDs(:,2)];
    as = [aAs(:,1), aBs(:,1), aCs(:,1), aDs(:,1),...
        aAs(:,2), aBs(:,2), aCs(:,2), aDs(:,2)];
    
    writetable(table(CD, 'VariableNames', {'CD'}), './outputs/length_cd.txt');
    writematrix(rs, './outputs/lib_pos.txt');
    writematrix(vs, './outputs/lib_vel.txt');
    writematrix(as, './outputs/lib_acc.txt');
    writematrix([omg2s(:,3), omg4s(:,3)], './outputs/lib_ang_vel.txt');
    writematrix([alp2s(:,3), alp4s(:,3)], './outputs/lib_ang_acc.txt');
    writematrix(Js, './outputs/lib_J.txt');
end