addpath('./functions/pva')
%Export XY for Pressure-S diagram of B
data = readtable('./outputs/Practical_diagram'); s = data{:,1}; P = 10*data{:,2};
% plot(s,P,'k')

HB = 39; alpha = 90*pi/180; beta = 90*pi/180; n1 = 500; iter = 20000;
AB = 104.75; 
CD = 85.653; HD = 39.0362; % iter = 20000
OA = HB/2; AC = OA;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;
dien_tich_piston=100; %cm^2
Q2=5*exp(1i*3*pi/2); %Newton
Q3=5*exp(1i*3*pi/2);
Q4=5*exp(1i*3*pi/2);
Q5=Q3;

so_bien=length(s);
phi = zeros(so_bien,1);
Mds=zeros(so_bien,1);

for count=1:so_bien
    rB = (OA+AB-s(count))*exp(1i*xOB);
    if (count>278)
        [~,~,xOA,~] = CPA4(OA, AB, rB);
        phi(count)=xOA;
    elseif (count>156 && count <=278)
        [xOA,~,~,~] = CPA4(OA, AB, rB);
        phi(count)=xOA;
    elseif (count>=36 && count<=156)
        [~,~,xOA,~] = CPA4(OA, AB, rB);
        phi(count)=xOA;
    elseif (count<36)
        [xOA,~,~,~] = CPA4(OA, AB, rB);
        phi(count)=xOA;
    end
end

phi_deg=phi*180/pi;
axis tight manual
% a=animate(phi,so_bien);       % test if phi_deg is true
    
for count=1:so_bien
    if (phi_deg(count)<0)
    phi_deg(count)=phi_deg(count)+360;
    end
end

theta=zeros(so_bien,1);
theta(1)=phi_deg(1);
for count=2:92
    theta(count)=theta(count-1)+abs(phi_deg(count)-phi_deg(count-1));
end

for count=93:140
    theta(count)=360+phi_deg(count);
end
for count=141:156
    theta(count)=theta(140);
end
for count =157:350
    theta(count)=theta(count-1)+abs(phi_deg(count)-phi_deg(count-1));
end
for count=351:356
    theta(count)=720+phi_deg(count);
end
F_B=P*dien_tich_piston/(10^4);
F_Bs=F_B*exp(1i*5*pi/4);
% plot(theta,F,'k')

i=1;
for count=1:so_bien
    rA = OA*exp(1i*phi(count));
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
    rD = OD*exp(1i*xOD);
    rcenterCA=(rC+rA)/2;
    rcenterAB=(rA+rB)/2;
    rS2=rcenterCA+(rcenterAB-rcenterCA)*AC/AB;
    rS4=(rC+rD)/2;
    
    vA = vel(rA, 0, omg1);
    [v3_T, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
    [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
    vB = vel(rB,v3_T, 0);
    vC = vel(rC, vC_T, omgC);
    [v5_T, omg4] = vCPA2(-vC, rD, rD-rC, 0, 0);
    vD = vel(rD, v5_T, 0);
    vS2=vel(rS2, 0, omg2);
    vS4=vel(rS4, 0, omg4);
    vB = vB/1000;
    vC = vC/1000;
    vD = vD/1000;
    vS2 = vS2/1000;
    vS4 = vS4/1000;
    one=cdot(F_Bs(i),vB);
    two=cdot(Q2,vS2);
    three=cdot(Q3,vB);
    four=cdot(Q4,vS2);
    five=cdot(Q5,vD);
    Mds(i)=(one+two+three+four+five)/omg1;
%     Mds(i)=abs(Mds(i));
    i=i+1;
    
end

% Ad=cumtrapz(Mds); %This starts from 45 degree xOA
%   Because Mds starts from 45 degree of xOA, so we will calculate
%   And plot from 0 to 720 degree 
theta_8pi = zeros(so_bien*2,1);
Md_8pi = zeros(so_bien*2,1);
for count=1:so_bien
    theta_8pi(count) = theta(count);
    theta_8pi(count+so_bien) = 720+theta(count);
    Md_8pi(count)=Mds(count);
    Md_8pi(count+so_bien)=Mds(count);
end
% plot(theta_8pi,Md_8pi);%This is right
theta_0_to_720 = zeros(so_bien,1);
Md_0_to_720 = zeros(so_bien,1);
for count=1:so_bien
    theta_0_to_720(count) = theta_8pi(count+92)-360;
    Md_0_to_720(count) = Md_8pi(count+92);
end
% plot(theta_0_to_720,Md_0_to_720,'k');

Ad=cumtrapz(Md_0_to_720);
% plot(theta_0_to_720,Ad);
Ad_energy = zeros(so_bien,2);
for count=1:so_bien
    Ad_energy(count,1) = theta_0_to_720(count);
    Ad_energy(count,2) = Ad(count);
end
% plot(theta,F_B);
% title ('F_B with respet to \phi_1');
% hold on
% figure
% plot(theta,Mds);
% title ('M_d with respet to \phi_1');
% figure
% plot(theta,Ad);
% title ('A_d with respet to \phi_1');

F_B_and_phi = zeros (so_bien,2);
for count=1:so_bien
    F_B_and_phi(count,1)=phi(count);
    F_B_and_phi(count,2)=F_B(count);
end
function P = cdot(F, v)
  P = real(F)*real(F) + imag(v)*imag(v);
end