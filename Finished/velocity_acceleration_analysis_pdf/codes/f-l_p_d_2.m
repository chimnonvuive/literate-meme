% find rE
eq1 = (rD(2)-y)/(rD(1)-x) == (rD(2)-rC(2))/(rD(1)-rC(1));
eq2 = (rC(1)-x)^2 + (rC(2)-y)^2 == CE^2;
sol = solve([eq1,eq2], [x,y]);
if subs(rC(1), theta, theta1) > subs(sol.x(1), theta, theta1)
    rE = subs(rE, [x, y], [sol.x(1), sol.y(1)]);
else
    rE = subs(rE, [x, y], [sol.x(2), sol.y(2)]);
end
end