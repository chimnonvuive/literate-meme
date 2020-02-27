function [vB, vC, vE, omg2, omg3] = d_vel(theta1)
global rD
[rB, rC, rE] = d_pos(theta1);
theta2 = atan((rC(2)-rB(2))/(rC(1)-rB(1)));
theta3 = atan((rC(2)-rD(2))/(rC(1)-rD(1)));
% find vB, vC, vE, omg2, omg3
vB = diff(rB); vC = diff(rC); vE = diff(rE);
omg2 = [0,0,diff(theta2)]; omg3 = [0,0,diff(theta3)];
end
function [aB, aC, aE, alp2, alp3] = d_acc(theta1)
[vB, vC, vE, omg2, omg3] = d_vel(theta1);
% find aB, aC, aE, alp2, alp3
aB = diff(vB); aC = diff(vC); aE = diff(vE);
alp2 = diff(omg2); alp3 = diff(omg3);
end