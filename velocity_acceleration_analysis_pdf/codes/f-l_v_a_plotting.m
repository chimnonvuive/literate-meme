global AB BC CD CE rD omg1 alp1
AB = .15; BC = .35; CD = .3; CE = .15;
rD = [.3,.3,0]; omg1 = [0,0,1]; alp1 = [0,0,0];

i = 0; phis = zeros(201, 1);
aCs = phis; vCs = phis;

for phi=0:(2*pi/200):2*pi
    i=i+1; phis(i)=phi;
    [~,vC,~,~,~] = vel(phi); vCs(i, 1) = vC(1);
    [~,aC,~,~,~] = acc(phi); aCs(i, 1) = aC(1);
end
plot(phis, vCs(:, 1))
plot(phis, aCs(:, 1))
