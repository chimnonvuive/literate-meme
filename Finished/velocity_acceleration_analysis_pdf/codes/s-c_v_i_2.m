% return the values
omg12 = double(subs(omg12, sol.x));
omg23 = double(subs(omg23, sol.y));
v23 = double(subs(v23, sol.z));
vC = -v23;
omg2 = omg1 + omg12;
end