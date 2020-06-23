function [aB, aC, aE, alp2, alp3] = c_acc(theta1)
global omg1 alp1 rD, syms x y z
% set variables
[rB, rC, rE] = pos(theta1);
[~,~,~, omg2, omg3] = vel(theta1);
aB = cross(alp1, rB) - omg1(3)^2*rB;
alp12 = [0,0,x]; alp23 = [0,0,y]; alp30 = [0,0,z];
% find alp12, alp23, alp30
eq1 = alp1 + alp12 + alp23 + alp30 == 0;
eq2 = cross(rB, alp12) + cross(rC, alp23) + cross(rD, alp30) - omg1(3)^2*rB - omg2(3)^2*(rC-rB) - omg3(3)^2*(rD-rC) == 0;
sol = solve([eq1(3),eq2(1),eq2(2)], [x,y,z]);