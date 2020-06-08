% Pressure acting on link 5
data = readtable('./outputs/phiD_PD.txt');
phis = data{:,1}; iter = length(phis);

rs = readtable('./outputs/lib_pos.txt');
rAs = rs{:,1} + rs{:,5}*1j;
rBs = rs{:,2} + rs{:,6}*1j;
rCs = rs{:,3} + rs{:,7}*1j;
rDs = rs{:,4} + rs{:,8}*1j;

vs = readtable('./outputs/lib_vel.txt');
vAs = vs{:,1} + vs{:,5}*1j;
vBs = vs{:,2} + vs{:,6}*1j;
vCs = vs{:,3} + vs{:,7}*1j;
vDs = vs{:,4} + vs{:,8}*1j;

as = readtable('./outputs/lib_acc.txt');
aAs = as{:,1} + as{:,5}*1j;
aBs = as{:,2} + as{:,6}*1j;
aCs = as{:,3} + as{:,7}*1j;
aDs = as{:,4} + as{:,8}*1j;

alps = readtable('./outputs/lib_ang_acc.txt');
alp2s = alps{:,1}; alp4s = alps{:,2};

% rA rotates from 135 deg to -45 deg
rA1 = [rAs(17501:20000);rAs(1:7500)];
rB1 = [rBs(17501:20000);rBs(1:7500)];
rC1 = [rCs(17501:20000);rCs(1:7500)];
rD1 = [rDs(17501:20000);rDs(1:7500)];

% rA rotates from -45 deg to 135 deg
rA2 = rAs(7501:17500);
rB2 = rBs(7501:17500);
rC2 = rCs(7501:17500);
rD2 = rDs(7501:17500);

% ----- Designed parameters -----
HB = 39; HD = 39.0074; OA = HB/2; AB = 104.395; AC = OA/2; CD = 94.8311;
alpha = 90; beta = 90; n1 = 500;
BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = deg2rad(90 - alpha/2); xOD = deg2rad(90 + alpha/2);
omg1 = convangvel(n1,'rpm','rad/s');

rAs_  = zeros(iter, 1); vAs_  = zeros(iter, 1); aAs_  = zeros(iter, 1);
rBs_  = zeros(iter, 1); vBs_  = zeros(iter, 1); aBs_  = zeros(iter, 1);
rCs_  = zeros(iter, 1); vCs_  = zeros(iter, 1); aCs_  = zeros(iter, 1);
rDs_  = zeros(iter, 1); vDs_  = zeros(iter, 1); aDs_  = zeros(iter, 1);
alp2s_ = zeros(iter, 1); alp4s_ = zeros(iter, 1);

% figure, sys(iter) = struct('cdata',[],'colormap',[]);

for i=1:iter
    dist = abs(rAs - OA*exp(1j*phis(i)));
    index = find(dist == min(dist));
    rAs_(i) = rAs(index); vAs_(i) = vAs(index); aAs_(i) = aAs(index);
    rBs_(i) = rBs(index); vBs_(i) = vBs(index); aBs_(i) = aBs(index);
    rCs_(i) = rCs(index); vCs_(i) = vCs(index); aCs_(i) = aCs(index);
    rDs_(i) = rDs(index); vDs_(i) = vDs(index); aDs_(i) = aDs(index);
    alp2s_(i) = alp2s(index); alp4s_(i) = alp4s(index);
%     plot_frame(rAs_(i), rBs_(i), rCs_(i), rDs_(i), sys(i));
end

rs_ = [rAs_, rBs_, rCs_, rDs_];
vs_ = [vAs_, vBs_, vCs_, vDs_];
as_ = [aAs_, aBs_, aCs_, aDs_];
writematrix([real(rs_), imag(rs_)], './outputs/eq_M_pos');
writematrix([real(vs_), imag(vs_)], './outputs/eq_M_vel');
writematrix([real(as_), imag(as_)], './outputs/eq_M_acc');
writematrix([alp2s_, alp4s_], './outputs/eq_M_alp');


function plot_frame(rA, rB, rC, rD, sys)

    link1 = [0, rA]; 
    link2_AB = [rA, rB];
    link2_AC = [rA, rC];
    link2_BC = [rB, rC]; 
    link4 = [rC, rD];
    
    hold on, axis equal
    plot(real(link1), imag(link1), 'bo-')
    plot(real(link2_AB), imag(link2_AB), 'ro-')
    plot(real(link2_AC), imag(link2_AC), 'ro-')
    plot(real(link2_BC), imag(link2_BC), 'ro-')
    plot(real(link4), imag(link4), 'go-')
    xlim([-150,150]); ylim([-150,150]);
    drawnow;
    sys = getframe;
    clf
end