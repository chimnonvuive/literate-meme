% find alp2, alp3
eq = aB + cross(alp2, rC-rB) - omg2(3)^2*(rC-rB) == aC;
sol = solve([eq(1),eq(2)], [x,y]);
% find aC, aE
aC = double(subs(aC, sol.y));
aE = double(subs(aE, sol.y));
alp2 = double(subs(alp2, sol.x));
alp3 = double(subs(alp3, sol.y));
end