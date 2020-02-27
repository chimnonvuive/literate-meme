function [vB, vC, omg2, v23] = c_vel(theta1)
global omg1, syms x y z
% set variables
[rB, rC] = pos(theta1); vB = cross(omg1, rB);
omg12 = [0,0,x]; omg23 = [0,0,y]; v23 = [z,0,0];
% find omg12, omg23, v23
eq1 = omg1 + omg12 + omg23 == 0;
eq2 = cross(rB, omg12) + cross(rC, omg23) + v23 == 0;
sol = solve([eq1(3),eq2(1),eq2(2)], [x,y,z]);