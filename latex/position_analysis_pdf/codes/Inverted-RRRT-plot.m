global AD BC AC rC
AD = 0.35; BC = 0.20; AC = 0.15; rC = [AC,0,0];
phi = pi/3; omg1 = [0,0,30*pi/30]; alp1 = [0,0,0];

hold on, axis equal
[rB, rD] = pos(pi/3);
plot([0, rD(1)], [0,rD(2)], 'r-o')
plot([rB(1), rC(1)], [rB(2), rC(2)], 'b-o');
