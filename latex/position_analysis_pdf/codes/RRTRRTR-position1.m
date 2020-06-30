function [rB, rD, rF, rG] = pos(phi1)

global AB DC DF AG rC

rB = AB*[cos(phi1), sin(phi1), 0];
pm = sign(rB(1));

% find rD
CB_angle = atan((rB(2)-rC(2)) / (rB(1)-rC(1)))+pi;
rD = double(rC + pm*DC*[cos(CB_angle), sin(CB_angle), 0]);



