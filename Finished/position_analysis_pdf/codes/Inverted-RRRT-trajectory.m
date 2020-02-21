global AD BC AC rC
AD = 0.35; BC = 0.20; AC = 0.15; rC = [AC,0,0];
phi = pi/3; omg1 = [0,0,30*pi/30]; alp1 = [0,0,0];

for phi=0:2*pi/32:2*pi
    hold on, axis equal
    [rB, rD] = pos(phi);
    plot([rC(1), rB(1)], [rC(2), rB(2)], 'b-o');
    plot([0, rD(1)], [0, rD(2)], 'r-o');
end
