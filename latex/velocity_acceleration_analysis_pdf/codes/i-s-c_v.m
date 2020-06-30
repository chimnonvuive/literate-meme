function [vB, vD, omg3, vB2B1] = vel(theta1)
global rC omg1, syms x y z
% set variables
[rB, rD] = pos(theta1);
vB2B1 = [x,y,0]; omg3 = [0,0,z];
vB = cross(omg3, rB-rC); vD = cross(omg1, rD);
% find vB2B1, omg3=omg2
eq1 = cross(omg3, rB-rC) == vB2B1 + cross(omg1, rB);
eq2 = dot(cross(cross(rB,[1,0,0]),rB), vB2B1) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);
% return the values
vB2B1 = double(subs(vB2B1, [x, y], [sol.x, sol.y]));
omg3 = double(subs(omg3, sol.z));
vB = double(subs(vB, sol.z));
end