function [aB, aC, aE, alp2, alp3] = acc(theta1)
global omg1 rD alp1, syms x y
% set variables
[rB, rC, rE] = pos(theta1);
[~,~,~, omg2, omg3] = vel(theta1);
alp2 = [0,0,x]; alp3 = [0,0,y];
aB = cross(alp1, rB) - omg1(3)^2*rB;
aC = cross(alp3, rC-rD) - omg3(3)^2*(rC-rD);
aE = cross(alp3, rE-rD) - omg3(3)^2*(rE-rD);