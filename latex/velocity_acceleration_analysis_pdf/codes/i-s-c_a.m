function [aB, aD, alp3] = acc(theta1)
global rC alp1 omg1, syms x y z
% set variables
[rB, rD] = pos(theta1); [~,~, omg3, vB2B1] = vel(theta1);
aB2B1 = [x,y,0]; alp3 = [0,0,z];
aB1 = cross(alp1, rB) - omg1(3)^2*rB;
aD = cross(alp1, rD) - omg1(3)^2*rD;
aB2 = cross(alp3, rB-rC) - omg3(3)^2*(rB-rC);