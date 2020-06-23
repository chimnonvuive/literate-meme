function [vB, vC, vE, omg2, omg3] = c_vel(theta1)
global rD omg1, syms x y z
% set variables
[rB, rC, rE] = pos(theta1); vB = cross(omg1, rB);
omg12 = [0,0,x]; omg23 = [0,0,y]; omg30 = [0,0,z];
% find omg12, omg23, omg30
eq1 = omg1 + omg12 + omg23 + omg30 == 0;
eq2 = cross(rB, omg12) + cross(rC, omg23) + cross(rD, omg30) == 0;
sol = solve([eq1(3),eq2(1),eq2(2)], [x,y,z]);