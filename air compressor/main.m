addpath('./functions'), addpath('./outputs')

HB = 39; alpha = 90*pi/180; beta = 90*pi/180; n1 = 500; iter = 20000;
AB = 104.75; 
CD = 85.653; HD = 39.0362; % iter = 20000
OA = HB/2; AC = OA;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;
dien_tich_piston=100; %cm^2
m2=0.5;
m3=m2; 
m4=m2; 
m5=m2;
Q2=m2*exp(1i*3*pi/2);
Q3=m3*exp(1i*3*pi/2);
Q4=m4*exp(1i*3*pi/2);
Q5=Q3;
J1 = 0.12;
J2 = J1; 
J3 = J1; 
J4 = J1;
% m1=6;
% m2=4;
% m3=1.5; 
% m4=2.5; 
% m5=1.5;
% Q2=m2*exp(1i*3*pi/2);
% Q3=m3*exp(1i*3*pi/2);
% Q4=m4*exp(1i*3*pi/2); 
% J1 = 0.5;
% J2 = 0.6; 
% J3 = 0.7; 
% J4 =0.8; 
delta=1/80;
muyE=1;
muyJ=1;

% For F_B, Ad
P_B = readtable('Practical diagram.csv'); s = P_B{:,1}; P = 10*P_B{:,2};
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
% a=animate(phi,so_bien);
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
F_B_and_phi = zeros (so_bien,2);
for count=1:so_bien
    F_B_and_phi(count,1)=phi(count);
    F_B_and_phi(count,2)=F_B(count);
end

%For J
ABj=AB/1000; HBj=HB/1000; CDj=CD/1000; HDj=HD/1000;
OAj=OA/1000;ACj=AC/1000;BCj=BC/1000;
Js = zeros(iter,1);
i=1;
for theta=linspace(0,2*pi,iter+1)
    rA = OAj*exp(1i*theta);
    [~,~,~,OB] = CPA2(xOB, ABj, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(ACj, BCj, rB-rA);
    rC = rA + ACj*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CDj, rC);
    rD = OD*exp(1i*xOD);
    rS2=(rC+rA)/2+(rB-rA)*AB/AC;
    rS4=(rC+rD)/2;
    

    vA = vel(rA, 0, omg1);
    [v3_T, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
    [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
    vC = vel(rC, vC_T, omgC);
    [v5_T, omg4] = vCPA2(-vC, rD, rD-rC, 0, 0);
    vS2=vel(rS2, 0, omg2);
    vS4=vel(rS4, 0, omg4);
    v3 = vel(rB, v3_T, 0);
    v5 = vel(rD, v5_T, 0);
    vS2s=abs(vS2);
    vS4s=abs(vS4);
    vBs = abs(v3);
    vDs = abs(v5);
    Js(i)=J1+J2*(omg2/omg1)^2+m2*(vS2s/omg1)^2+m3*(vBs/omg1)^2+...
        J4*(omg4/omg1)^2+m4*(vS4s/omg1)^2+m5*(vDs/omg1)^2;
    i = i + 1;
end
Js_0_to_720 = zeros(iter*2+2,2);
for count=1:iter+1
    Js_0_to_720(count,1)=(count-1)*360/(iter);
    Js_0_to_720(count+iter+1,1)=360+(count-1)*360/(iter);
    Js_0_to_720(count,2)=Js(count);
    Js_0_to_720(count+iter+1,2)=Js(count);
end
% hold on
% plot(linspace(0,360,iter+1),Js,'k');
% plot(linspace(360,720,iter+1),Js,'k');

% For F_D, -Ac
data = readtable('P at D diagram.csv'); s_D = data{:,1}; P_D = data{:,2};
load('rAs_.mat', 'rAs_');
load('rBs_.mat', 'rBs_');
load('rCs_.mat', 'rCs_');
load('rDs_.mat', 'rDs_');
so_bien2=length(rAs_);
phi = zeros(so_bien2,1);
Mcs=zeros(so_bien2,1);
for count=1:so_bien2
    phi(count)=angle(rAs_(count));
end
% Algorithm to reverse array
phi_reverse = zeros(so_bien2,1);
for count=1:163
    phi_reverse(count) = phi(326+1-count);
    phi_reverse(326+1-count) = phi(count);
end
phi_deg=phi_reverse*180/pi;
for count=1:so_bien2
    if (phi_deg(count)<0)
    phi_deg(count)=phi_deg(count)+360;
    end
end
theta=zeros(so_bien2,1);
theta(1)=phi_deg(1);
for count=2:209
    theta(count)=theta(count-1)+abs(phi_deg(count)-phi_deg(count-1));
end
theta=zeros(so_bien2,1);
theta(1)=phi_deg(1);
for count=2:209
    theta(count)=theta(count-1)+abs(phi_deg(count)-phi_deg(count-1));
end
for count=210:so_bien2
    theta(count)=360+phi_deg(count);
end
F_D=P_D*dien_tich_piston/(10^4);
F_Ds=F_D*exp(1i*5*pi/4);
i=1;
for count=1:so_bien2
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
theta8pi = zeros(so_bien2*4,1);
Mc8pi = zeros(so_bien2*4,1);
for count=1:so_bien2
    theta8pi(count) = theta(count);
    theta8pi(count+so_bien2) = 360+theta(count);
    theta8pi(count+so_bien2*2) = 720+theta(count);
    theta8pi(count+so_bien2*3) = 1080+theta(count);
    Mc8pi(count) = Mcs(count);
    Mc8pi(count+so_bien2) = Mcs(count);
    Mc8pi(count+so_bien2*2) = Mcs(count);
    Mc8pi(count+so_bien2*3) = Mcs(count);
end
% plot(theta4pi,Mcs4pi);
theta_0_to_7202 = zeros(so_bien2*2,1);
Mc_0_to_7202 = zeros(so_bien2*2,1);
for count=1:so_bien2*2
    theta_0_to_7202(count) = theta8pi(count+535)-720;
    Mc_0_to_720(count) = Mc8pi(count+535);
end
% plot(theta_0_to_7202,Mc_0_to_720,'k');
Ac=cumtrapz(Mc_0_to_720);
% plot(theta_0_to_7202,Ac);
%h1 = length of last point in Ac diagram
%h2 = length of last point in Ad diagram
% TLX_P_D_real= h2/h1*TLX_P_D_fake;
h1=Ac(length(Ac));
h2=Ad(length(Ad));
P_D_real=h2*P_D/h1;
F_D_real=P_D_real*dien_tich_piston/(10^4);
F_D_real_0to720=zeros(so_bien2*2,1);
Mc_0_to_720_real=Mc_0_to_720*(h2/h1);
Ac_real=cumtrapz(Mc_0_to_720_real);
% plot(theta_0_to_7202,Ac_real);
Ac_energy = zeros(so_bien2*2,2);
for count=1:2*so_bien2
    Ac_energy(count,1) = theta_0_to_7202(count);
    Ac_energy(count,2) = Ac_real(count);
end
F_D_and_phi = zeros (so_bien2*2,2);
for count=1:so_bien2
    F_D_and_phi(count,1)=phi(count);
    F_D_and_phi(count+so_bien2,1)=phi(count);
    F_D_and_phi(count,2)=F_D(count);
    F_D_and_phi(count+so_bien2,2)=F_D(count);
end

%For energy calculation
% plot(Ac_energy(:,1),Ac_energy(:,2),'k',Ad_energy(:,1),Ad_energy(:,2),'b');
new_Ad= zeros (length(Ac_energy),1);
[x, index] = unique(Ad_energy(:,1)); 
for count=1:length(Ac_energy)
    new_Ad(count)=interp1(x,Ad_energy(index,2),Ac_energy(count,1),'pchip');
end
% plot(Ac_energy(:,1),new_Ad,'o-k',Ac_energy(:,1),Ac_energy(:,2),'o-b');
E=new_Ad-Ac_energy(:,2)+1/2*omg1^2*Js_0_to_720(1,2);
% plot(Ac_energy(:,1),E,'k');
new_Js= zeros (length(E),1);
[x, index] = unique(Js_0_to_720(:,1)); 
for count=1:length(E)
    new_Js(count) = interp1(x,Js_0_to_720(index,2),Ac_energy(count,1),'pchip');
end
hold on
plot(new_Js,E,'k');
plot([new_Js(1) new_Js(652)],[E(1) E(652)],'k');
%chon delta(he so khong deu cho phep)=1/80
%Chon muy J =1; muy E=1;
psi_max=atan(muyJ/(2*muyE))*omg1^2*(1+delta/2)^2;
psi_min=atan(muyJ/(2*muyE))*omg1^2*(1-delta/2)^2;
%Find the right tangent
%y=ax+b with a=tan(psi) and b is the point on E axis 
% (a is psi_maxor and b is psi_min)
pointA=zeros(length(E),1);
pointB=zeros(length(E),1);
for count=1:length(E)
    pointA(count)=E(count)-tan(psi_max)*new_Js(count);
    pointB(count)=E(count)-tan(psi_max)*new_Js(count);
end 
a=max(pointA);
b=min(pointB);
%Plot to show 2 tangent line
x_td = zeros (26,1);
for i=1:26
    x_td(i)=-0.1+(i-1)*0.01;
end
ytmax=tan(psi_max)*x_td+a;
ytmin=tan(psi_min)*x_td+b;
plot(x_td,ytmax,'r',x_td,ytmin,'r');
grid on
Jd=muyJ*(ytmax(11)-ytmin(11))/(tan(psi_max)-tan(psi_min))-Js_0_to_720(1,2);