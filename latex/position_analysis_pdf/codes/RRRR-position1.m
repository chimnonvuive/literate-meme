function [rB, rC, rE] = pos(phi1)
global AB BC CD CE rD
syms x y

rB = AB*[cos(phi1), sin(phi1), 0];

% find rC
rC = [x, y, 0];
eq1 = (x-rB(1))^2 + (y-rB(2))^2 == BC^2;
eq2 = (x-rD(1))^2 + (y-rD(2))^2 == CD^2;
sol = solve(eq1, eq2, x, y);
if 0>sol.x(1)
    rC = double(subs(rC, [x, y], [sol.x(1), sol.y(1)]));
else
    rC = double(subs(rC, [x, y], [sol.x(2), sol.y(2)])); end
