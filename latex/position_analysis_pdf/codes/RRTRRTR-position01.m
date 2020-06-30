function [rB, rD, rF, rG] = pos(phi1)

global AB DC DF AG rC
syms x y

rB = AB*[cos(phi1), sin(phi1), 0];

% find rD
rD = [x, y, 0];
eq1 = (x-rC(1))^2 + (y - rC(2))^2 == DC^2;
eq2 = (y-rC(2))/(x-rC(1)) == (y-rB(2))/(x-rB(1));
sol = solve(eq1, eq2, x, y);

