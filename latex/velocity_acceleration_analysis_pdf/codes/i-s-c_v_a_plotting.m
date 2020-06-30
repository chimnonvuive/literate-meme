global AD BC AC rC
AD = 0.35; BC = 0.20; AC = 0.15; rC = [AC,0,0];
theta1 = pi/3; omg1 = [0,0,pi]; alp1 = [0,0,0];

i=0; theta1s = zeros(201, 1);
vBs = theta1s; aBs = theta1s;

for theta1=0:(2*pi/200):2*pi
    i=i+1; theta1s(i) = theta1;
    vB2 = vel(theta1); vBs(i) = vB2(1);
    aB2 = acc(theta1); aBs(i) = aB2(1);
end
plot(theta1s, vBs)
plot(theta1s, aBs)