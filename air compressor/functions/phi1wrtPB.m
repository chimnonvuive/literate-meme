function [rAs_, rBs_, rCs_, rDs_, vAs_, vBs_, vCs_, vDs_, aAs_, aBs_, aCs_, aDs_, alp2s_, alp4s_] = phi1wrtPB(sB)
    
%     input: sB
%     output: corresponding system position rAs_, rBs_, rCs_, rDs_
% 
%     minimum rB at period 1 and 2, which means max sB at period 1 and 2
%     the strategy is as follows:
%     find the middle point "a" of data sB as a relative point separating 2
%     periods approximately. After separation, determine the maximum
%     distance sB of 2 periods. Use the 2 points to narrow down the trough
%     between them, which is min_mid
    a = fix(length(sB)/2); period1 = sB(1:a); period2 = sB(a+1:end);
    crest_1 = find(period1 == max(period1),1);
    crest_2 = find(period2 == max(period2),1) + a;

    % find trough of sB between 2 crests
    trough = sB(crest_1+1:crest_2);
    trough = find(trough == min(trough),1) + crest_1;

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
    
    % specify data point at singularity
    max_lens = length(rAs); xOB = angle(rBs(1));
    s1 = fix(xOB/(2*pi)*max_lens);
    s2 = fix((xOB+pi)/(2*pi)*max_lens);
    
    % rA rotates from 45 deg to -135 deg ccw
    % from 0 to T/2 and T to 3T/2, we use this segment
    rA1 = rAs(s1+1:s2);
    rB1 = rBs(s1+1:s2);
    rC1 = rCs(s1+1:s2);
    rD1 = rDs(s1+1:s2);
    
    vA1 = vAs(s1+1:s2);
    vB1 = vBs(s1+1:s2);
    vC1 = vCs(s1+1:s2);
    vD1 = vDs(s1+1:s2);
    
    aA1 = aAs(s1+1:s2);
    aB1 = aBs(s1+1:s2);
    aC1 = aCs(s1+1:s2);
    aD1 = aDs(s1+1:s2);
    
    alp21 = alp2s(s1+1:s2);
    alp41 = alp4s(s1+1:s2);
    
    % rA rotates from -135 deg to 45 deg ccw
    % from T/2 to T and 3T/2 to 2T, we use this segment
    rA2 = [rAs(s2+1:max_lens);rAs(1:s1)];
    rB2 = [rBs(s2+1:max_lens);rBs(1:s1)];
    rC2 = [rCs(s2+1:max_lens);rCs(1:s1)];
    rD2 = [rDs(s2+1:max_lens);rDs(1:s1)];
    
    vA2 = [vAs(s2+1:max_lens);vAs(1:s1)];
    vB2 = [vBs(s2+1:max_lens);vBs(1:s1)];
    vC2 = [vCs(s2+1:max_lens);vCs(1:s1)];
    vD2 = [vDs(s2+1:max_lens);vDs(1:s1)];
    
    aA2 = [aAs(s2+1:max_lens);aAs(1:s1)];
    aB2 = [aBs(s2+1:max_lens);aBs(1:s1)];
    aC2 = [aCs(s2+1:max_lens);aCs(1:s1)];
    aD2 = [aDs(s2+1:max_lens);aDs(1:s1)];
    
    alp22 = [alp2s(s2+1:max_lens);alp2s(1:s1)];
    alp42 = [alp4s(s2+1:max_lens);alp4s(1:s1)];
    
    % ----- Designed parameters -----
    OA = abs(rAs(1)); AB = abs(rAs(1)-rBs(1)); 

    iter = length(sB);
    rAs_ = zeros(iter, 1);
    rBs_ = zeros(iter, 1);
    rCs_ = zeros(iter, 1);
    rDs_ = zeros(iter, 1);
    vAs_ = zeros(iter, 1);
    vBs_ = zeros(iter, 1);
    vCs_ = zeros(iter, 1);
    vDs_ = zeros(iter, 1);
    aAs_ = zeros(iter, 1);
    aBs_ = zeros(iter, 1);
    aCs_ = zeros(iter, 1);
    aDs_ = zeros(iter, 1);
    alp2s_ = zeros(iter, 1);
    alp4s_ = zeros(iter, 1);

    for i=1:iter
        if (1<=i && i<crest_1) || (trough<=i && i<crest_2)
            dist = abs(rB1 - (OA+AB-sB(i))*exp(1j*xOB));
            index = find(dist == min(dist));
            rAs_(i) = rA1(index); rBs_(i) = rB1(index);
            rCs_(i) = rC1(index); rDs_(i) = rD1(index);
            vAs_(i) = vA1(index); vBs_(i) = vB1(index);
            vCs_(i) = vC1(index); vDs_(i) = vD1(index);
            aAs_(i) = aA1(index); aBs_(i) = aB1(index);
            aCs_(i) = aC1(index); aDs_(i) = aD1(index);
            alp2s_(i) = alp21(index); alp4s_(i) = alp41(index);
        else
            dist = abs(rB2 - (OA+AB-sB(i))*exp(1j*xOB));
            index = find(dist == min(dist));
            rAs_(i) = rA2(index); rBs_(i) = rB2(index);
            rCs_(i) = rC2(index); rDs_(i) = rD2(index);
            vAs_(i) = vA2(index); vBs_(i) = vB2(index);
            vCs_(i) = vC2(index); vDs_(i) = vD2(index);
            aAs_(i) = aA2(index); aBs_(i) = aB2(index);
            aCs_(i) = aC2(index); aDs_(i) = aD2(index);
            alp2s_(i) = alp22(index); alp4s_(i) = alp42(index);
        end
    end
end
