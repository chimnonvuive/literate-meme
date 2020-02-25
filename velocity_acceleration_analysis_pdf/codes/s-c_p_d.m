function [rB, rC] = d_pos(phi)
global AB BC, syms x theta(t)

rB = [AB*cos(theta(t)), AB*sin(theta(t)), 0];
rC = [x,0,0];

% find rC
eq = (rB(1)-x)^2 + rB(2)^2 == BC^2;
sol = solve(eq, x);
if subs(sol(1), theta, phi) > subs(rB(1), theta, phi)
    rC = subs(rC, x, sol(1));
else
    rC = subs(rC, x, sol(2)); end
end
