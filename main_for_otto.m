alpha = 90*pi/180; beta = 90*pi/180; n1 = 500; iter = 20000;
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;
HB = 39; AB = 104.75; CD = 85.653; HD = 39.0362;
OA = HB/2; AC = OA;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
dien_tich_piston=100; %cm^2
m2=0.5; Q2=m2*exp(1i*3*pi/2);   J1 = 0.12;
m3=m2;  Q3=m3*exp(1i*3*pi/2);   J2 = J1; 
m4=m2;  Q4=m4*exp(1i*3*pi/2);   J3 = J1; 
m5=m2;  Q5=Q3;                  J4 = J1;
delta=1/80; muyE=1; muyJ=1;

% For F_B, Ad
P_B = readtable('practical_otto_diagram.csv'); s_B = P_B{:,1}; PofB = P_B{:,2};
so_bien=length(s_B);
phi1 = zeros(so_bien,1);
Mds=zeros(so_bien,1);
for count=1:so_bien
    rB = (OA+AB-s_B(count))*exp(1i*xOB);
    if (count>=600 && count <=so_bien)
        [~,~,xOA,~] = CPA4(OA, AB, rB);
        phi1(count)=xOA;
    elseif (count>=380 && count <=599)
        [xOA,~,~,~] = CPA4(OA, AB, rB);
        phi1(count)=xOA;
    elseif (count>=171 && count <=379)
        [~,~,xOA,~] = CPA4(OA, AB, rB);
        phi1(count)=xOA;
    elseif (count>=1 && count <=170)
        [xOA,~,~,~] = CPA4(OA, AB, rB);
        phi1(count)=xOA;
    end    
end 
deg1=double(phi1)*180/pi;
% a=animate(phi,so_bien);
for count=1:so_bien
    if (deg1(count)<0)
    deg1(count)=deg1(count)+360;
    end
end
for count=317:723       
    deg1(count)=deg1(count)+360; end
for count=724:so_bien   
    deg1(count)=deg1(count)+720; end
F_B=PofB*dien_tich_piston/(10^4);
FB_8pi=zeros(so_bien*2,1);
for count=1:so_bien
    FB_8pi(count)=F_B(count);
    FB_8pi(count+so_bien)=F_B(count);
end
F_Bs=F_B*exp(1i*5*pi/4);
i=1;
for count=1:so_bien
    rA = OA*exp(1i*phi1(count));
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
    rD = OD*exp(1i*xOD);
    centerCA=(rC+rA)/2;
    centerAB=rA*0.37/OA; %Using Solidworks to file center of mass
    rS2=centerCA+(centerAB-centerCA)*AC/AB;
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
    Mds(i)=(cdot(F_Bs(i),vB)+cdot(Q2,vS2)+cdot(Q2,vS2)+cdot(Q4,vS2)+cdot(Q5,vD))/omg1;
    i=i+1;
end
deg1_8pi = zeros(so_bien*2,1);
Md_8pi = zeros(so_bien*2,1);
for count=1:so_bien
    deg1_8pi(count) = deg1(count);
    deg1_8pi(count+so_bien) = 720+deg1(count);
    Md_8pi(count)=Mds(count);
    Md_8pi(count+so_bien)=Mds(count);
end
deg1_unique=zeros(720,1);
for count=1:720
    deg1_unique(count)=count;
end
Md_unique= zeros (length(deg1_unique),1);
[deg1_fix, index] = unique(deg1_8pi); 
for count=1:length(deg1_unique)
    Md_unique(count)=interp1(deg1_fix,Md_8pi(index),deg1_unique(count),'pchip');
end
FB_unique= zeros (length(deg1_unique),1);
[deg1_fix2, index] = unique(deg1_8pi-720); 
for count=1:length(deg1_unique)
    FB_unique(count)=interp1(deg1_fix2,FB_8pi(index),deg1_unique(count),'pchip');
end

Ad=cumtrapz(Md_unique);
Ad_deg1 = zeros(length(deg1_unique),2);
for count=1:length(deg1_unique)
    Ad_deg1(count,1) = deg1_unique(count);
    Ad_deg1(count,2) = Ad(count);
end
FB_deg1 = zeros (length(deg1_unique),2);
for count=1:length(deg1_unique)
    FB_deg1(count,1)=deg1_unique(count);
    FB_deg1(count,2)=FB_deg1(count);
end

%For J
Js = zeros(iter,1);
i=1;
for theta=linspace(0,2*pi,iter+1)
    rA = (OA/1000)*exp(1i*theta);
    [~,~,~,OB] = CPA2(xOB, AB/1000, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC/1000, BC/1000, rB-rA);
    rC = rA + (AC/1000)*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD/1000, rC);
    rD = OD*exp(1i*xOD);
    centerCA=(rC+rA)/2;
    centerAB=rA*0/OA; %Using Solidworks to file center of mass
    rS2=centerCA+(centerAB-centerCA)*AC/AB;
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
Js4pi = zeros(iter*2+2,2);
for count=1:iter+1
    Js4pi(count,1)=(count-1)*360/(iter);
    Js4pi(count+iter+1,1)=360+(count-1)*360/(iter);
    Js4pi(count,2)=Js(count);
    Js4pi(count+iter+1,2)=Js(count);
end

% For F_D, -Ac
P_D = readtable('P at D diagram.csv'); s_D = P_D{:,1}; PofD = P_D{:,2};
load('rAs_.mat', 'rAs_');
load('rBs_.mat', 'rBs_');
load('rCs_.mat', 'rCs_');
load('rDs_.mat', 'rDs_');
so_bien2=length(rAs_);
phi2 = zeros(so_bien2,1);
Mcs=zeros(so_bien2,1);
for count=1:so_bien2
    phi2(count)=angle(rAs_(count));
end
deg2=phi2*180/pi;
for count=262:so_bien2
    deg2(count)=-180+deg2(count)-180;
end
for count=1:so_bien2
    deg2(count)=135+135-deg2(count);
end
F_D=PofD*dien_tich_piston/(10^4);
F_Ds=F_D*exp(1i*5*pi/4);
i=1;
for count=1:so_bien2
    rA = OA*exp(1i*phi2(count));
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
deg2_8pi = zeros(so_bien2*4,1);
Mc_8pi = zeros(so_bien2*4,1);
for count=1:so_bien2
    deg2_8pi(count) = deg2(count);
    deg2_8pi(count+so_bien2) = 360+deg2(count);
    deg2_8pi(count+so_bien2*2) = 720+deg2(count);
    deg2_8pi(count+so_bien2*3) = 1080+deg2(count);
    Mc_8pi(count) = Mcs(count);
    Mc_8pi(count+so_bien2) = Mcs(count);
    Mc_8pi(count+so_bien2*2) = Mcs(count);
    Mc_8pi(count+so_bien2*3) = Mcs(count);
end
deg2_8pi=deg2_8pi-360;
Mc_unique= zeros (length(deg1_unique),1);
[deg2_fix, index] = unique(deg2_8pi); 
for count=1:length(deg1_unique)
    Mc_unique(count)=interp1(deg2_fix,Mc_8pi(index),deg1_unique(count),'pchip');
end
Ac=cumtrapz(Mc_unique);
% plot(deg1_unique,Ac);
%h1 = length of last point in Ac diagram
%h2 = length of last point in Ad diagram
% TLX_P_D_real= h2/h1*TLX_P_D_fake;
h1=Ac(length(Ac));
h2=Ad(length(Ad));
PD_real=h2*PofD/h1;
FD_real=PD_real*dien_tich_piston/(10^4);
FD_real_8pi=zeros(so_bien2*4,1);
for count=1:so_bien2
    FD_real_8pi(count) = FD_real(count);
    FD_real_8pi(count+so_bien2) = FD_real(count);
    FD_real_8pi(count+so_bien2*2) = FD_real(count);
    FD_real_8pi(count+so_bien2*3) = FD_real(count);
end
FD_unique= zeros (length(deg1_unique),1);
for count=1:length(deg1_unique)
    FD_unique(count)=interp1(deg2_fix,FD_real_8pi(index),deg1_unique(count),'pchip');
end
Mc_real=Mc_unique*(h2/h1);
Ac_real=cumtrapz(Mc_real);
Ac_deg1 = zeros(length(deg1_unique),2);
for count=1:length(deg1_unique)
    Ac_deg1(count,1) = deg1_unique(count);
    Ac_deg1(count,2) = Ac_real(count);
end
FD_deg1 = zeros (length(deg1_unique),2);
for count=1:length(deg1_unique)
    FD_deg1(count,1)=deg1_unique(count);
    FD_deg1(count,2)=FD_unique(count);
end

%For energy calculation
E=Ad_deg1(:,2)-Ac_deg1(:,2)+1/2*omg1^2*Js4pi(1,2);
Js_unique= zeros (length(deg1_unique),1);
[deg1_Js, index] = unique(Js4pi(:,1)); 
for count=1:length(deg1_unique)
    Js_unique(count) = interp1(deg1_Js,Js4pi(index,2),deg1_unique(count),'pchip');
end
%chon delta(he so khong deu cho phep)=1/80
% Chon muy J =1; muy E=1;
psi_max=atan(muyJ/(2*muyE))*omg1^2*(1+delta/2)^2;
psi_min=atan(muyJ/(2*muyE))*omg1^2*(1-delta/2)^2;
%Find the right tangent
%y=ax+b with a=tan(psi) and b is the point on E axis 
% (a is psi_maxor and b is psi_min)
pointA=zeros(length(E),1);
pointB=zeros(length(E),1);
for count=1:length(E)
    pointA(count)=E(count)-tan(psi_max)*Js_unique(count);
    pointB(count)=E(count)-tan(psi_max)*Js_unique(count);
end 
a=max(pointA);
b=min(pointB);
x_td = zeros (26,1);
for i=1:26
    x_td(i)=-0.1+(i-1)*0.01;
end
ytmax=tan(psi_max)*x_td+a;
ytmin=tan(psi_min)*x_td+b;
Jd=muyJ*(ytmax(11)-ytmin(11))/(tan(psi_max)-tan(psi_min))-Js4pi(1,2);

% FOR PB_s and PD_s
% plot(s_B,PofB);
% plot(s_D,PofD);

% FOR FB and FD
% plot(deg1_unique,FB_unique);
% plot(deg1_unique,FD_unique);

% FOR MOMENT DONG, MOMENT CAN
% plot(deg1_unique,Md_unique);
% plot(deg1_unique,Mc_unique);

% FOR CONG DONG, CONG CAN
% plot(Ad_deg1(:,1),Ad_deg1(:,2),'b',Ac_deg1(:,1),Ac_deg1(:,2),'r');

% FOR E-J and flywheel
% plot(deg1_unique,E,'k');
% grid on, hold on
% plot(Js_unique,E,'k');
% plot(x_td,ytmax,'r',x_td,ytmin,'r');

function P = cdot(F, v)
  P = real(F)*real(F) + imag(v)*imag(v);
end