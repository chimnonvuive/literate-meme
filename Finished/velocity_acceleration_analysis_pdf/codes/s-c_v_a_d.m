function [vB, vC, omg2] = d_vel(theta1)
[rB, rC] = d_pos(theta1); theta2 = atan((rC(2)-rB(2))/(rC(1)-rB(1)));
% find vB, vC, omg2
vB = diff(rB); vC = diff(rC); omg2 = [0,0,diff(theta2)];
end
function [aB, aC, alp2] = d_acc(theta1)
[vB, vC, omg2] = d_vel(theta1);
% find aB, aC, alp2
aB = diff(vB); aC = diff(vC); alp2 = diff(omg2);
end