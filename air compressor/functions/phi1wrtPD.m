function phi1 = phi1wrtPD(sD, OA, AB, alpha)

    xOD = pi/2 + alpha/2; 

    crest = find(sD == max(sD),1,'first');

    rs = readtable('outputs/lib_pos.txt');
    rAs = rs{:,1} + rs{:,5}*1j;
    rBs = rs{:,2} + rs{:,6}*1j;
    rCs = rs{:,3} + rs{:,7}*1j;
    rDs = rs{:,4} + rs{:,8}*1j;
    
    len_rs = length(rAs);
    dp1 = fix(xOD*len_rs/2/pi);
    dp2 = fix((xOD+pi)*len_rs/2/pi);
    % rA rotates from 135 deg to -45 deg
    rA1 = [rAs(dp2+1:end);rAs(1:dp1)];
    rB1 = [rBs(dp2+1:end);rBs(1:dp1)];
    rC1 = [rCs(dp2+1:end);rCs(1:dp1)];
    rD1 = [rDs(dp2+1:end);rDs(1:dp1)];
    
    % rA rotates from -45 deg to 135 deg
    rA2 = rAs(dp1+1:dp2);
    rB2 = rBs(dp1+1:dp2);
    rC2 = rCs(dp1+1:dp2);
    rD2 = rDs(dp1+1:dp2);

    iter = length(sD);
    rAs_ = zeros(iter, 1);
    rBs_ = zeros(iter, 1);
    rCs_ = zeros(iter, 1);
    rDs_ = zeros(iter, 1);

    for i=1:iter
        if i > crest
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
    phi1 = angle(rAs_);
%     animation(rAs_,rBs_,rCs_,rDs_);
end
