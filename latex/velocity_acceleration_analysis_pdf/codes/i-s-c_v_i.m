function [vB, vD, omg3, v12] = c_vel(theta1)
global rC omg1, syms x y z t
% set variables
[rB, rD] = pos(theta1);
omg23 = [0,0,x]; omg30 = [0,0,y]; v12 = [z,t,0];
% find omg23, omg30, a12
eq1 = omg1 + omg23 + omg30 == 0;
eq2 = cross(rB, omg23) + cross(rC, omg30) + v12 == 0;
eq3 = dot(cross(cross(rB, [1,0,0]), rB), v12) == 0;
sol = solve([eq1(3), eq2(1), eq2(2), eq3], [x,y,z,t]);
