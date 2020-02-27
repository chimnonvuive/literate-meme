global AB BC omg1 alp1
AB = .5; BC = 1; omg1 = [0,0,1]; alp1 = [0,0,-1];

i = 0; theta1s = zeros(201, 1);
vCs = theta1s; aCs = theta1s;

for theta1=0:(2*pi/200):2*pi
    i=i+1; theta1s(i)=phi;
    [~,vC] = vel(theta1); vCs(i) = vC(1);
    [~,aC] = acc(theta1); aCs(i) = aC(1);
end
plot(theta1s, vCs)
plot(theta1s, aCs)
