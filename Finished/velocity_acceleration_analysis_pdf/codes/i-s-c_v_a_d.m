function [vB, vD, omg3] = d_vel(phi)
global rC
[rB, rD] = d_pos(phi);
% find vB, vD, omg3
vB = diff(rB); vD = diff(rD);
phi3 = atan((rB(2)-rC(2))/(rB(1)-rC(1)));
omg3 = [0, 0, diff(phi3)];
end


function [aB, aD, alp3] = d_acc(phi)
[vB, vD, omg3] = d_vel(phi);
% find aB, aD, alp3
aB = diff(vB); aD = diff(vD); alp3 = diff(omg3);
end