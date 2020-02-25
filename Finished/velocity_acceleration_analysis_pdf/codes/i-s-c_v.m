function [vB2, vD, omg3, vB2B1] = vel(phi)
global rC omg1, syms x y z

% set variables
[rB, rD] = pos(phi); vB2B1 = [x,y,0]; omg3 = [0,0,z];
vB2 = cross(omg3, rB-rC); vB1 = cross(omg1, rB);

% find vB2B1, omg3=omg2 and return the values
eq1 = vB2 == vB2B1 + vB1;
eq2 = dot(cross(cross(rB, [1,0,0]), rB), vB2B1) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);

vB2B1 = double(subs(vB2B1, [x, y], [sol.x, sol.y]));
omg3 = double(subs(omg3, z, sol.z));
vB2 = double(subs(vB, z, sol.z)); vD = cross(omg1, rD); end