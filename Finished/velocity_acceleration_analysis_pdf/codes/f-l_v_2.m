% find omg2, omg3
eq = vB + cross(omg2, rC-rB) == vC;
sol = solve([eq(1),eq(2)], [x,y]);
% return the values
omg2 = double(subs(omg2, sol.x));
omg3 = double(subs(omg3, sol.y));
vC = double(subs(vC, sol.y));
vE = double(subs(vE, sol.y));
end