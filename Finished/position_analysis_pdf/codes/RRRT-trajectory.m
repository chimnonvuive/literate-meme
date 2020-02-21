global AB BC omg1 alp1
AB = .5; BC = 1; omg1 = [0,0,1]; alp1 = [0,0,-1];

for phi=0:2*pi/32:2*pi
    hold on, axis equal
    [rB, rC] = pos(phi);
    plot([0, rB(1)], [0,rB(2)], 'r-o');
    plot([rB(1), rC(1)], [rB(2), rC(2)], 'b-o');
end