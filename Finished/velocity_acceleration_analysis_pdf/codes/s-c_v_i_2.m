
% return the values
omg12 = double(subs(omg12, x, sol.x));
omg23 = double(subs(omg23, y, sol.y));
v23 = double(subs(v23, z, sol.z));
vB = cross(omg1, rB);
vC = -v23; omg2 = omg1 + omg12;
end