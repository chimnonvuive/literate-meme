
% return the values
omg23 = double(subs(omg23, x, sol.x));
omg30 = double(subs(omg30, y, sol.y));
v12 = double(subs(v12, [z,t], [sol.z,sol.t]));
omg3 = omg1 + omg23;
vB2 = cross(omg3, rB-rC);
vD = cross(omg1, rD);
end