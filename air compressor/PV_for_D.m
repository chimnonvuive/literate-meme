%PV diagram for piston D
data = readtable('P at D diagram.csv'); s_D = data{:,1}; P_D = data{:,2};
load('rAs_.mat', 'rAs_');
load('rBs_.mat', 'rBs_');
load('rCs_.mat', 'rCs_');
load('rDs_.mat', 'rDs_');
% plot(s_D,P_D)

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

so_bien=length(rAs_);
phi = zeros(so_bien,1);
Mcs=zeros(so_bien,1);

for count=1:so_bien
    phi(count)=angle(rAs_(count));
end

% Algorithm to reverse array
phi_reverse = zeros(326,1);
for count=1:163
    phi_reverse(count) = phi(326+1-count);
    phi_reverse(326+1-count) = phi(count);
end

phi_deg=phi_reverse*180/pi;
for count=1:so_bien
    if (phi_deg(count)<0)
    phi_deg(count)=phi_deg(count)+360;
    end
end

theta=zeros(so_bien,1);
theta(1)=phi_deg(1);
for count=2:209
    theta(count)=theta(count-1)+abs(phi_deg(count)-phi_deg(count-1));
end

for count=210:326
    theta(count)=360+phi_deg(count);
end

F_D=P_D*dien_tich_piston/(10^4);
F_Ds=F_D*exp(1i*5*pi/4);

% plot(theta,F_D,'k');
% hold on
% plot(theta+360,F_D,'k');

% a=animate(phi,so_bien);

% for count=1:so_bien
%     hold on, axis equal
%     axis tight manual
%     plot([0 cos(phi_reverse(count))],[0 sin(phi_reverse(count))],'k');
%     text(cos(phi_reverse(count)),sin(phi_reverse(count)), [num2str(count)]);
% end

i=1;
for count=1:so_bien
    rA = OA*exp(1i*phi(count));
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
    rD = OD*exp(1i*xOD);
    
    vA = vel(rA, 0, omg1);
    [v3_T, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
    [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
    vB = vel(rB,v3_T, 0);
    vC = vel(rC, vC_T, omgC);
    [v5_T, omg4] = vCPA2(-vC, rD, rD-rC, 0, 0);
    vD = vel(rD, v5_T, 0);
    vD = vD/1000;
    Mcs(i)=cdot(F_Ds(i),vD)/omg1;
    i=i+1;
end
%Theta only 360 degree => calculate theta2 and Mcs to have 720 degree
theta8pi = zeros(so_bien*4,1);
Mc8pi = zeros(so_bien*4,1);
for count=1:so_bien
    theta8pi(count) = theta(count);
    theta8pi(count+so_bien) = 360+theta(count);
    theta8pi(count+so_bien*2) = 720+theta(count);
    theta8pi(count+so_bien*3) = 1080+theta(count);
    Mc8pi(count) = Mcs(count);
    Mc8pi(count+so_bien) = Mcs(count);
    Mc8pi(count+so_bien*2) = Mcs(count);
    Mc8pi(count+so_bien*3) = Mcs(count);
end
% plot(theta4pi,Mcs4pi);

theta_0_to_720 = zeros(so_bien*2,1);
Mc_0_to_720 = zeros(so_bien*2,1);
for count=1:so_bien*2
    theta_0_to_720(count) = theta8pi(count+535)-720;
    Mc_0_to_720(count) = Mc8pi(count+535);
end
% plot(theta_0_to_720,Mc_0_to_720,'k');
Ac=cumtrapz(Mc_0_to_720);
% plot(theta_0_to_720,Ac);

%h2 = length of last point in Ac diagram
%h1 = length of last point in Ad diagram
% TLX_P_D_real= h2/h1*TLX_P_D_fake;
h1=4.568799755737407e+08;
h2=4.746806066538444e+06;
P_D_real=h2*P_D/h1;
F_D_real=P_D_real*dien_tich_piston/(10^4);
max_Md_0_to_720 = 5.051579358628707e+04;
max_Mc_0_to_720 = 2.148591741247515e+06;
Mc_0_to_720_real=Mc_0_to_720*(h2/h1);
Ac_real=cumtrapz(Mc_0_to_720_real);
plot(theta_0_to_720,Ac_real);

Ac_energy = zeros(so_bien*2,2);
for count=1:2*so_bien
    Ac_energy(count,1) = theta_0_to_720(count);
    Ac_energy(count,2) = Ac_real(count);
end

F_D_and_phi = zeros (so_bien*2,2);
for count=1:so_bien
    F_D_and_phi(count,1)=phi(count);
    F_D_and_phi(count+so_bien,1)=phi(count);
    F_D_and_phi(count,2)=F_D(count);
    F_D_and_phi(count+so_bien,2)=F_D(count);
end

function P = cdot(F, v)
  P = real(F)*real(F) + imag(v)*imag(v);
end