function [vB, vC, omg2] = vel(theta1)
global omg1, syms x y
% set variables
[rB, rC] = pos(theta1);
vB = cross(omg1, rB);
vC = [x,0,0]; omg2 = [0,0,y];
% find vC, omg2
eq = vB + cross(omg2, rC-rB) == vC;
sol = solve([eq(1),eq(2)], [x,y]);
% return the values
vC = double(subs(vC, sol.x));
omg2 = double(subs(omg2, sol.y));
end