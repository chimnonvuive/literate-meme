% find rE
rE = [x, y, 0];
eq1 = (rD(2)-rE(2))/(rD(1)-rE(1)) == (rD(2)-rC(2))/(rD(1)-rC(1));
eq2 = (rC(1)-rE(1))^2 + (rC(2)-rE(2))^2 == CE^2;
sol = solve(eq1, eq2, x, y);
if rC(1)>double(sol.x(1))
    rE = double(subs(rE, [x, y], [sol.x(1), sol.y(1)]));
else
    rE = double(subs(rE, [x, y], [sol.x(2), sol.y(2)])); end
end
