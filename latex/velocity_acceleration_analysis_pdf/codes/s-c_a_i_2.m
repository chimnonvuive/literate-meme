% return the values
alp12 = double(subs(alp12, sol.x));
alp23 = double(subs(alp23, sol.y));
a23 = double(subs(a23, sol.z));
aC = -a23;
alp2 = alp1 + alp12;
end