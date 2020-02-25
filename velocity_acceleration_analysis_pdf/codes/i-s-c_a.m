function [aB2, aD, alp2] = acc(phi)
global rC alp1 omg1, syms x y z

% set variables
[rB, rD] = pos(phi); [~,~, omg3, vB2B1] = vel(phi);
aB2B1r = [x,y,0]; alp2 = [0,0,z];
aB1 = cross(alp1, rB) - omg1(3)^2*rB;
aB2 = cross(alp2, rB-rC) - omg3(3)^2*(rB-rC);

% find aB2B1r, alp3=alp2
eq1 = aB2 == aB1 + aB2B1r + 2*cross(omg1, vB2B1);
eq2 = dot(cross(cross(rB, [1,0,0]), rB), aB2B1r) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);