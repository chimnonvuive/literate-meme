function [rB, rD] = d_pos(theta1)
global AD BC rC, syms l_AB theta(t)
% set variables
rB = l_AB*[cos(theta(t)), sin(theta(t)), 0];
rD = AD*[cos(theta(t)), sin(theta(t)), 0];
% find rB
eq = (rB(1)-rC(1))^2 + (rB(2)-rC(2))^2 == BC^2;
sol = solve(eq, l_AB);
if double(subs(sol(1), theta, theta1))>0
    rB = subs(rB, l_AB, sol(1));
else
    rB = subs(rB, l_AB, sol(2)); end
end