% return the values
alp12 = double(subs(alp12, sol.x));
alp23 = double(subs(alp23, sol.y));
alp30 = double(subs(alp30, sol.z));
alp2 = alp1 + alp12; alp3 = -alp30;
aC = cross(alp3, rC-rD) - omg3(3)^2*(rC-rD);
aE = cross(alp3, rE-rD) - omg3(3)^2*(rE-rD);
end