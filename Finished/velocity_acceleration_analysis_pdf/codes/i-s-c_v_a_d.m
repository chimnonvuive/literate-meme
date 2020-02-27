function [vB, vD, omg3] = d_vel(theta1)
global rC
[rB, rD] = d_pos(theta1);
theta3 = atan((rB(2)-rC(2))/(rB(1)-rC(1)));
% find vB, vD, omg3
vB = diff(rB); vD = diff(rD);
omg3 = [0,0,diff(theta3)];
end
function [aB, aD, alp3] = d_acc(theta1)
[vB, vD, omg3] = d_vel(theta1);
% find aB, aD, alp3
aB = diff(vB); aD = diff(vD); alp3 = diff(omg3);
end