rD1 = double(subs(rD, [sol.x(1), sol.y(1)]));
rD2 = double(subs(rD, [sol.x(2), sol.y(2)]));
if rD1(1) <= rC(1), rD = rD1; else, rD = rD2; end

% find rF, rG
rFC_angle = atan((rB(2)-rC(2))/(rB(1)-rC(1)));
rG_angle = atan(rD(2)/rD(1)) + pi;
rF = double(rD + DF*[cos(rFC_angle), sin(rFC_angle), 0]);
rG = double(AG*[cos(rG_angle), sin(rG_angle), 0]);

end
