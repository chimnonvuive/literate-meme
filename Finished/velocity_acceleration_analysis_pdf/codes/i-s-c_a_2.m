
% return the values
aB2B1r = double(subs(aB2B1, [x,y], [sol.x, sol.y]));
aB2 = double(subs(aB2, sol.z));
aD = cross(alp1, rD) - omg1(3)^2*rD;
alp2 = double(subs(alp2, sol.z));
end