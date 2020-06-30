function [rB, rC] = d_pos(theta1)
global AB BC, syms x theta(t)
% set variables
rB = AB*[cos(theta(t)), sin(theta(t)), 0];
rC = [x,0,0];
% find rC
eq = (rB(1)-x)^2 + rB(2)^2 == BC^2;
sol = solve(eq, x);
if subs(sol(1),theta,theta1) > subs(sol(2),theta,theta1)
    rC = subs(rC, sol(1));
else
    rC = subs(rC, sol(2));
end
end