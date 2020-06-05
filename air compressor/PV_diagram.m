%Export XY for Pressure-S diagram
TB = readtable('s_p.csv'); sB = TB{:,1}; PB = TB{:,2};
% plot(s,P,'k')

HB = 39; alpha = 90*pi/180; beta = 90*pi/180; n1 = 500; %iter = 20000;
AB = 104.75; 
CD = 85.653; HD = 39.0362; % iter = 20000
OA = HB/2; AC = OA;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;
S_piston = 100; %cm^2
Q2 = 5*exp(1i*3*pi/2); Q3 = Q2; Q4 = Q2; Q5 = Q2; 

iter = length(sB);
% rAs = zeros(iter,1);
% rBs = zeros(iter,1);
% rCs = zeros(iter,1); 
% rDs = zeros(iter,1);
vBs  = zeros(iter,1);
vCs  = zeros(iter,1); 
vDs  = zeros(iter,1);
vS2s = zeros(iter,1);
vS4s = zeros(iter,1);
Mds  = zeros(iter,1);
phis = zeros(iter,1);
thetas = zeros(iter,1);

for i=1:iter
    rB = (OA+AB-sB(i))*exp(1i*xOB);
    [xOA,~,~,~] = CPA4(OA, AB, rB);
    phis(i)=xOA;
end

thetas(1) = phis(1);

for i=1:iter-1
    thetas(i+1)=thetas(i)+abs(phis(i+1)-phis(i));
end

F_B = PB*S_piston/10000;
F_Bs = F_B*exp(1i*pi/4);
% plot(theta,F,'k')

for i=1:iter
    rA = OA*exp(1i*phis(i));
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
    rD = OD*exp(1i*xOD);
    
    vA = vel(rA, 0, omg1);
    [v3_T, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
    [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
    vBs(i) = vel(rB,v3_T, 0);
    vCs(i) = vel(rC, vC_T, omgC);
    [v5_T, omg4] = vCPA2(-vCs(i), rD, rD-rC, 0, 0);
    vDs(i) = vel(rD, v5_T, 0);
    vS2s(i) = (vA+vBs(i))/2;
    vS4s(i) = (vCs(i)+vDs(i))/2;
    vBs(i) = vBs(i)/1000;
    vCs(i) = vCs(i)/1000;
    vDs(i) = vDs(i)/1000;
    vS2s(i) = vS2s(i)/1000;
    vS4s(i) = vS4s(i)/1000;
    Mds(i)=(cdot(F_Bs(i),vBs(i)) + cdot(Q2,vS2s(i)) + cdot(Q3,vBs(i)) + cdot(Q4,vS4s(i)) + cdot(Q5,vDs(i)))/omg1;
end