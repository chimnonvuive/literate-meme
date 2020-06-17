function [phi1] = phi1wrtPD(sD)

    crest = find(sD == max(sD),1,'first');

    rs = readtable('outputs/lib_pos.txt');
    rAs = rs{:,1} + rs{:,5}*1j;
    rBs = rs{:,2} + rs{:,6}*1j;
    rCs = rs{:,3} + rs{:,7}*1j;
    rDs = rs{:,4} + rs{:,8}*1j;

%     vs = readtable('outputs/lib_vel.txt');
%     vAs = vs{:,1} + vs{:,5}*1j;
%     vBs = vs{:,2} + vs{:,6}*1j;
%     vCs = vs{:,3} + vs{:,7}*1j;
%     vDs = vs{:,4} + vs{:,8}*1j;
% 
%     as = readtable('outputs/lib_acc.txt');
%     aAs = as{:,1} + as{:,5}*1j;
%     aBs = as{:,2} + as{:,6}*1j;
%     aCs = as{:,3} + as{:,7}*1j;
%     aDs = as{:,4} + as{:,8}*1j;
% 
%     alps = readtable('outputs/lib_ang_acc.txt');
%     alp2s = alps{:,1}; alp4s = alps{:,2};


    % rA rotates from 135 deg to -45 deg
%     rA1 = [rAs(17501:20000);rAs(1:7500)]/1e3;
%     rB1 = [rBs(17501:20000);rBs(1:7500)]/1e3;
%     rC1 = [rCs(17501:20000);rCs(1:7500)]/1e3;
%     rD1 = [rDs(17501:20000);rDs(1:7500)]/1e3;
%     
%     rA1 = [rAs(12501:20000);rAs(1:2500)]/1e3;
%     rB1 = [rBs(12501:20000);rBs(1:2500)]/1e3;
%     rC1 = [rCs(12501:20000);rCs(1:2500)]/1e3;
%     rD1 = [rDs(12501:20000);rDs(1:2500)]/1e3;
% 
%     % rA rotates from -45 deg to 135 deg
%     rA2 = rAs(2501:12500)/1e3;
%     rB2 = rBs(2501:12500)/1e3;
%     rC2 = rCs(2501:12500)/1e3;
%     rD2 = rDs(2501:12500)/1e3;
    
    len_rs = length(rAs);
    xOD = angle(rDs(1));
    dp1 = fix(xOD*len_rs/2/pi);
    dp2 = fix((xOD+pi)*len_rs/2/pi);
    % rA rotates from 135 deg to -45 deg
    rA1 = [rAs(dp2+1:end);rAs(1:dp1)]/1e3;
    rB1 = [rBs(dp2+1:end);rBs(1:dp1)]/1e3;
    rC1 = [rCs(dp2+1:end);rCs(1:dp1)]/1e3;
    rD1 = [rDs(dp2+1:end);rDs(1:dp1)]/1e3;
    
    % rA rotates from -45 deg to 135 deg
    rA2 = rAs(dp1+1:dp2)/1e3;
    rB2 = rBs(dp1+1:dp2)/1e3;
    rC2 = rCs(dp1+1:dp2)/1e3;
    rD2 = rDs(dp1+1:dp2)/1e3;
    
    

    % ----- Designed parameters -----
    OA = abs(rAs(1))/1e3; AB = abs(rAs(1)-rBs(1))/1e3;

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
