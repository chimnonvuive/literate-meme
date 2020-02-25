global AD BC AC rC
AD = 0.35; BC = 0.20; AC = 0.15; rC = [AC,0,0];
phi = pi/3; omg1 = [0,0,pi]; alp1 = [0,0,0];

i = 0; phis = zeros(201, 1);
vBs = phis; aBs = phis;

for phi=0:(2*pi/200):2*pi
    i=i+1; phis(i)=phi;
    [vB2,~,~,~] = vel(phi); vBs(i, 1) = vB2(1);
    [aB2,~,~]   = acc(phi); aBs(i, 1) = aB2(1);
end
plot(phis, vBs(:, 1))
plot(phis, aBs(:, 1))