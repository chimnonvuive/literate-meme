function [aB, aC, alp2] = c_acc(phi)
global omg1 alp1, syms x y z

% set variables
[rB, rC] = pos(phi); [~,~,omg2,~] = c_vel(phi);
alp12 = [0,0,x]; alp23 = [0,0,y]; a23 = [z,0,0];

% find aC, alp2
eq1 = alp1 + alp12 + alp23 == 0;
eq2 = cross(rB, alp12) + cross(rC, alp23) + a23 - omg1(3)^2*rB - omg2(3)^2*(rC-rB) == 0;
sol = solve([eq1(3), eq2(1), eq2(2)], [x, y, z]);