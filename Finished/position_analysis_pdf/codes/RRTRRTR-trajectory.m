global AB AC DC DF AG rA rC omg1 alp1
AB = .15; AC = .1; DC = .15; DF = .4; AG = .3;
rA = [0,0,0]; rC = [0,AC,0];
omg1 = [0,0,5/3*pi]; alp1 = [0,0,0];

for phi=0:2*pi/32:2*pi
    hold on, axis equal
    [rB, rD, rF, rG] = pos(phi);
    plot([rA(1), rB(1)], [rA(2), rB(2)], 'r-o');
    plot([rD(1), rF(1)], [rD(2), rF(2)], 'g-o');
    plot([rA(1), rG(1)], [rA(2), rG(2)], 'b-o');
end
