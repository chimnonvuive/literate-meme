function [vB, vC, vE, omg2, omg3] = vel(theta1)
global omg1 rD, syms x y
% set variables
[rB, rC, rE] = pos(theta1);
omg2 = [0,0,x]; omg3 = [0,0,y];
vB = cross(omg1, rB);
vC = cross(omg3, rC-rD);
vE = cross(omg3, rE-rD);