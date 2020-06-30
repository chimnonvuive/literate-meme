global AB BC CD CE rD omg1 alp1
AB = .15; BC = .35; CD = .3; CE = .15;
rD = [.3 .3 0]; omg1 = [0, 0, 2*pi]; alp1 = [0,0,0];

hold on
[rB, rC, rE] = pos(pi/4);
plot([0, rB(1)], [0,rB(2)], 'r-o')
plot([rB(1), rC(1)], [rB(2), rC(2)], 'g-o')
plot([rD(1), rE(1)], [rD(2), rE(2)], 'b-o')
axis equal
