% Pressure acting on link 5
TD = readtable('B_pos.csv'); sD = TD{:,1};
index_max_sD = find(sD == max(sD));

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

iter = length(sD);
rAs_ = zeros(iter, 1);
rBs_ = zeros(iter, 1);
rCs_ = zeros(iter, 1);
rDs_ = zeros(iter, 1);
phis = zeros(iter, 1);
for i=1:iter
    if i > index_max_sD
        dist = abs(rD1 - (OA+AB-sD(i))*exp(1j*xOD));
        index = find(dist == min(dist));
        phis(i) = angle(rA1(index));
        rAs_(i) = rA1(index); rBs_(i) = rB1(index);
        rCs_(i) = rC1(index); rDs_(i) = rD1(index);
    else
        dist = abs(rD2 - (OA+AB-sD(i))*exp(1j*xOD));
        index = find(dist == min(dist));
        phis(i) = angle(rA2(index));
        rAs_(i) = rA2(index); rBs_(i) = rB2(index);
        rCs_(i) = rC2(index); rDs_(i) = rD2(index);
    end
end

animation(rAs_,rBs_,rCs_,rDs_)
