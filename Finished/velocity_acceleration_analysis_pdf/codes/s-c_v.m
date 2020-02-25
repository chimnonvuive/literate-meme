function [vB, vC, omg2] = vel(phi)
global omg1, syms x y

% set variables
[rB, rC] = pos(phi); vB = cross(omg1, rB);
vC = [x,0,0]; omg2 = [0,0,y];

% find vC, omg2
eq = vB + cross(omg2, rC-rB) == vC;
sol = solve([eq(1), eq(2)], [x, y]);

% return the values
vC = double(subs(vC, x, sol.x));
omg2 = double(subs(omg2, y, sol.y));
end
