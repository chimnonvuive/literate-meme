function [aB, aD, alp3] = c_acc(theta1)
global rC omg1 alp1, syms x y z t
% set variables
[rB, rD] = pos(theta1); [vB, vD, omg2, v12] = c_vel(theta1);
alp23 = [0,0,x]; alp30 = [0,0,y]; a12 = [z,t,0];
% find alp23, alp30, a12
eq1 = alp1 + alp23 + alp30 == 0;
eq2 = cross(rB, alp23) + cross(rC, alp30) + a12 + 2*cross(omg1, v12) - omg2(3)^2*(rC-rB) - omg1(3)^2*rB == 0;
eq3 = dot(cross(cross(rB, [1,0,0]), rB), a12) == 0;
sol = solve([eq1(3), eq2(1), eq2(2), eq3], [x,y,z,t]);