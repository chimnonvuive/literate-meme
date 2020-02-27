% return the results
alp23 = double(subs(alp23, x, sol.x));
alp30 = double(subs(alp30, y, sol.y));
a12 = double(subs(a12, [z,t], [sol.z,sol.t]));
alp3 = alp1 + alp23;
aB = cross(alp3, rB-rC)-omg2(3)^2*(rB-rC);
aD = cross(alp1, rD)-omg1(3)^2*rD;
end