function [aB, aC, alp2] = acc(theta1)
global omg1 alp1, syms x y
% set variables
[rB, rC] = pos(theta1);
[~,~, omg2] = vel(theta1);
aB = cross(alp1, rB) - norm(omg1)^2*rB;
aC = [x,0,0]; alp2 = [0,0,y];
% find aC, alp2
eq = aB + cross(alp2, rC-rB) - omg2(3)^2*(rC-rB) == aC;
sol = solve(eq(1), eq(2), x, y);
% return the values
aC = double(subs(aC, sol.x));
alp2 = double(subs(alp2, sol.y));
end