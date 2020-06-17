function [rAs_, rBs_, rCs_, rDs_] = phi1wrtPD(sD)

    index_max_sD = find(sD == max(sD));

    rs = readtable('outputs/lib_pos.txt');
    rAs = rs{:,1} + rs{:,5}*1j;
    rBs = rs{:,2} + rs{:,6}*1j;
    rCs = rs{:,3} + rs{:,7}*1j;
    rDs = rs{:,4} + rs{:,8}*1j;

    vs = readtable('outputs/lib_vel.txt');
    vAs = vs{:,1} + vs{:,5}*1j;
    vBs = vs{:,2} + vs{:,6}*1j;
    vCs = vs{:,3} + vs{:,7}*1j;
    vDs = vs{:,4} + vs{:,8}*1j;

    as = readtable('outputs/lib_acc.txt');
    aAs = as{:,1} + as{:,5}*1j;
    aBs = as{:,2} + as{:,6}*1j;
    aCs = as{:,3} + as{:,7}*1j;
    aDs = as{:,4} + as{:,8}*1j;

    alps = readtable('outputs/lib_ang_acc.txt');
    alp2s = alps{:,1}; alp4s = alps{:,2};


    % rA rotates from 135 deg to -45 deg
    rA1 = [rAs(17501:20000);rAs(1:7500)];
    rB1 = [rBs(17501:20000);rBs(1:7500)];
    rC1 = [rCs(17501:20000);rCs(1:7500)];
    rD1 = [rDs(17501:20000);rDs(1:7500)];

    % rA rotates from -45 deg to 135 deg
    rA2 = rAs(7501:17500);
    rB2 = rBs(7501:17500);
    rC2 = rCs(7501:17500);
    rD2 = rDs(7501:17500);

    % ----- Designed parameters -----
    OA = abs(rAs(1)); AB = abs(rAs(1)-rBs(1)); xOD = angle(rDs(1));

    iter = length(sD);
    rAs_ = zeros(iter, 1);
    rBs_ = zeros(iter, 1);
    rCs_ = zeros(iter, 1);
    rDs_ = zeros(iter, 1);

    for i=1:iter
        if i > index_max_sD
            dist = abs(rD1 - (OA+AB-sD(i))*exp(1j*xOD));
            index = find(dist == min(dist));
            rAs_(i) = rA1(index); rBs_(i) = rB1(index);
            rCs_(i) = rC1(index); rDs_(i) = rD1(index);
        else
            dist = abs(rD2 - (OA+AB-sD(i))*exp(1j*xOD));
            index = find(dist == min(dist));
            rAs_(i) = rA2(index); rBs_(i) = rB2(index);
            rCs_(i) = rC2(index); rDs_(i) = rD2(index);
        end
    end
end
