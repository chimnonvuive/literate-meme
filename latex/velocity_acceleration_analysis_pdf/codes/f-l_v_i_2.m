% return the values
omg12 = double(subs(omg12, sol.x));
omg23 = double(subs(omg23, sol.y));
omg30 = double(subs(omg30, sol.z));
omg2 = omg1 + omg12; omg3 = -omg30;
vC = cross(omg3, rC-rD);
vE = cross(omg3, rE-rD);
end