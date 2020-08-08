function [phi1] = phi1wrtPB(sB, OA, AB, alpha)
    
    rs = readtable('outputs/lib_pos.txt');
    
    xOB = pi/2 - alpha/2;

    mid_sB = fix(length(sB)/2);
    a = sB(1:mid_sB); b = sB(mid_sB+1:end);
    crest1 = find(a == max(a),1,'first');
    crest2 = find(b == max(b),1,'first') + mid_sB;
    c = sB(crest1:crest2);
    trough = find(c == min(c),1,'first') + crest1;
    
    rAs = rs{:,1} + rs{:,5}*1j;
    rBs = rs{:,2} + rs{:,6}*1j;
    rCs = rs{:,3} + rs{:,7}*1j;
    rDs = rs{:,4} + rs{:,8}*1j;

    len_rs = length(rAs);
    dp1 = fix(xOB*len_rs/2/pi);
    dp2 = fix((xOB+pi)*len_rs/2/pi);
    
    % rA rotates from 45 deg to -135 deg
    rA1 = rAs(dp1+1:dp2); length(rA1);
    rB1 = rBs(dp1+1:dp2);
    rC1 = rCs(dp1+1:dp2);
    rD1 = rDs(dp1+1:dp2);
    
    % rA rotates from -135 deg to 45 deg
    rA2 = [rAs(dp2+1:end);rAs(1:dp1)];
    rB2 = [rBs(dp2+1:end);rBs(1:dp1)];
    rC2 = [rCs(dp2+1:end);rCs(1:dp1)];
    rD2 = [rDs(dp2+1:end);rDs(1:dp1)];

    iter = length(sB);
    rAs_ = zeros(iter, 1);
    rBs_ = zeros(iter, 1);
    rCs_ = zeros(iter, 1);
    rDs_ = zeros(iter, 1);

    for i=1:iter
        if 0<i && i<=crest1 || trough<i && i<=crest2
            dist = abs(rB1 - (OA+AB-sB(i))*exp(1j*xOB));
            index = find(dist == min(dist));
            rAs_(i) = rA1(index); rBs_(i) = rB1(index);
            rCs_(i) = rC1(index); rDs_(i) = rD1(index);
        else
            dist = abs(rB2 - (OA+AB-sB(i))*exp(1j*xOB));
            index = find(dist == min(dist));
            rAs_(i) = rA2(index); rBs_(i) = rB2(index);
            rCs_(i) = rC2(index); rDs_(i) = rD2(index);
        end
    end
    phi1 = angle(rAs_);
end
