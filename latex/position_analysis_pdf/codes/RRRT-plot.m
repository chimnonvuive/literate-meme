global AB BC omg1 alp1
syms theta(t)

AB = .5; BC = 1; omg1 = [0,0,1]; alp1 = [0,0,-1];

hold on
[rB, rC] = pos(pi/3);
plot([0, rB(1)], [0,rB(2)], 'r-o');
plot([rB(1), rC(1)], [rB(2), rC(2)], 'b-o');
axis equal