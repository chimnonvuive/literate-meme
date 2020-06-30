% find aB2B1, alp3=alp2
eq1 = aB2 == aB1 + aB2B1 + 2*cross(omg1, vB2B1);
eq2 = dot(cross(cross(rB,[1,0,0]),rB), aB2B1) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);
% return the values
aB2B1 = double(subs(aB2B1, [x,y], [sol.x, sol.y]));
aB = double(subs(aB2, sol.z));
alp3 = double(subs(alp3, sol.z));
end