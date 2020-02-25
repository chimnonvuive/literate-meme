
% return the values
alp12 = double(subs(alp12, x, sol.x));
alp23 = double(subs(alp23, y, sol.y));
a23 = double(subs(a23, z, sol.z));
aB = cross(alp1, rB) - omg1(3)^2*rB;
aC = -a23; alp2 = alp1 + alp12;
end