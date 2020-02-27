global AB BC CD CE rD omg1 alp1
AB = .15; BC = .35; CD = .3; CE = .15;
rD = [.3,.3,0]; omg1 = [0,0,1]; alp1 = [0,0,0];

i=0; theta1s = zeros(201,1);
aCs = theta1s; vCs = theta1s;

for theta1=0:(2*pi/200):2*pi
    i=i+1;
    [~,vC] = vel(theta1); vCs(i) = vC(1);
    [~,aC] = acc(theta1); aCs(i) = aC(1);
    theta1s(i)=theta1;
end
plot(theta1s, vCs)
plot(theta1s, aCs)