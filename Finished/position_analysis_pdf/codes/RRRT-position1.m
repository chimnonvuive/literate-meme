function [rB, rC] = pos(phi1)
global AB BC
syms x

rB = AB*[cos(phi1), sin(phi1), 0]; rC = [x, 0, 0];

% find rC
eq = (rB(1)-rC(1))^2 + rB(2)^2 == BC^2;
sol = solve(eq, x);
rC1 = subs(rC, sol(1)); rC2 = subs(rC, sol(2));
if rC1(1) > rB(1), rC = double(rC1);
else, rC = double(rC2); end
end
