function [rB, rD] = pos(phi1)
global AD BC rC
syms l_AB

rB = l_AB*[cos(phi1), sin(phi1), 0];
rD = AD*[cos(phi1), sin(phi1), 0];

% find rB
eq = (rB(1)-rC(1))^2 + (rB(2)-rC(2))^2 == BC^2;
sol = solve(eq, l_AB);
if sol(1)>0, rB = double(subs(rB, sol(1)));
else, rB = double(subs(rB, sol(2))); end
end
