alpha = 90*pi/180; beta = 60*pi/180; n1 = 1700; iter = 20000;
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;
HB = 39e-3; AB = 104.75e-3; CD = 85.653e-3; HD = 39.0362e-3;
OA = HB/2; AC = OA;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
A_piston=100e-4; %m^2
g=9.81;
m2=4;    Q2=m2*g*exp(1i*3*pi/2);   J1 = 0.07;
m3=1.5;  Q3=m3*g*exp(1i*3*pi/2);   J2 = 0.12;
m4=2.5;  Q4=m4*g*exp(1i*3*pi/2); 
m5=1.5;  Q5=m3*g*exp(1i*3*pi/2);   J4 = 0.08;
delta=1/80; muyE=1; muyJ=1;


% For F_B, Ad
P_B = readtable('practical_otto_diagram.csv');
s_B = P_B{:,1}/1000; PofB = P_B{:,2};
num=length(s_B);
vBnum=zeros(num,2);
vDnum=zeros(num,2);
phi1 = zeros(num,1);
Mds=zeros(num,1);
for i=1:num
    rB = (OA+AB-s_B(i))*exp(1i*xOB);
    if (i>=600 && i <=num)
        [~,~,xOA,~] = CPA4(OA, AB, rB);
        phi1(i)=xOA;
    elseif (i>=380 && i <=599)
        [xOA,~,~,~] = CPA4(OA, AB, rB);
        phi1(i)=xOA;
    elseif (i>=171 && i <=379)
        [~,~,xOA,~] = CPA4(OA, AB, rB);
        phi1(i)=xOA;
    elseif (i>=1 && i <=170)
        [xOA,~,~,~] = CPA4(OA, AB, rB);
        phi1(i)=xOA;
    end
end 
deg1=double(real(phi1));
% a=animate(phi,num);
for i=1:num
    if (deg1(i)<0)
    deg1(i)=deg1(i)+2*pi;
    end
end
for i=317:723       
    deg1(i)=deg1(i)+2*pi; end
for i=724:num   
    deg1(i)=deg1(i)+4*pi; end
F_B=PofB*A_piston;
FB_8pi=[F_B;F_B];
F_Bs=F_B*exp(1i*5*pi/4);
for i=1:num
    rA = OA*exp(1i*phi1(i));
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
    rD = OD*exp(1i*xOD);
    centerCA=(rC+rA)/2;
    centerAB=(rA+rB)/2;
    rS2=centerCA+(centerAB-centerCA)*AC/AB;
    rS4=(rC+rD)/2;   
    vA = vel(rA, 0, omg1);
    [v3_T, omg2] = vCPA2(-vA, rB, rB-rA, 0, 0);
    [vC_T, omgC] = vCPA1(-vA, rC, rC-rA, 0, omg2);
    vB = vel(rB,v3_T, 0);   vBnum(i,1)=real(vB); vBnum(i,2)=imag(vB);
    vC = vel(rC, vC_T, omgC);
    [v5_T, omg4] = vCPA2(-vC, rD, rD-rC, 0, 0);
    vD = vel(rD, v5_T, 0); vDnum(i,1)=real(vD); vDnum(i,2)=imag(vD);
    vS2=vel(rS2, 0, omg2);
    vS4=vel(rS4, 0, omg4);
    Mds(i)=(cdot(F_Bs(i),vB)+cdot(Q3,vB)+cdot(Q2,vS2)+cdot(Q4,vS2)+cdot(Q5,vD))/omg1;
end
deg1_8pi = zeros(num*2,1);
Md_8pi=[Mds;Mds];
for i=1:num
    deg1_8pi(i) = real(deg1(i));
    deg1_8pi(i+num) = 4*pi+deg1(i);
end
deg1_unique=linspace(0,4*pi,720);
[deg1_fix, index] = unique(deg1_8pi-4*pi);
Md_unique=interp1(deg1_fix,Md_8pi(index),deg1_unique,'pchip');
[deg1_fix2, index] = unique(deg1_8pi-4*pi); 
FB_unique=interp1(deg1_fix2,FB_8pi(index),deg1_unique,'pchip');
Ad=cumtrapz(deg1_unique(2),Md_unique);
Ad_deg1 = zeros(length(deg1_unique),2);
for i=1:length(deg1_unique)
    Ad_deg1(i,1) = deg1_unique(i);
    Ad_deg1(i,2) = Ad(i);
end
FB_deg1 = zeros (length(deg1_unique),2);
for i=1:length(deg1_unique)
    FB_deg1(i,1)=deg1_unique(i);
    FB_deg1(i,2)=FB_unique(i);
end

% %For J
Js = zeros(iter,1); theta = 0;
for i=1:iter+1
    theta = 2*pi/(iter+1)*i;
    rA = (OA)*exp(1i*theta);
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + (AC)*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
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
end
Js4pi = zeros(iter*2+2,2);
for i=1:iter+1
    Js4pi(i,1)=(i-1)*360/(iter);
    Js4pi(i+iter+1,1)=360+(i-1)*360/(iter);
    Js4pi(i,2)=Js(i);
    Js4pi(i+iter+1,2)=Js(i);
end

% % For F_D, -Ac
P_D = readtable('P at D diagram.csv'); s_D = P_D{:,1}/1000; PofD = P_D{:,2};
load('rAs_.mat', 'rAs_');
load('rBs_.mat', 'rBs_');
load('rCs_.mat', 'rCs_');
load('rDs_.mat', 'rDs_');
num2=length(rAs_);
phi2 = zeros(num2,1);
Mcs=zeros(num2,1);
for i=1:num2
    phi2(i)=angle(rAs_(i));
end
for i=262:num2
    phi2(i)=-pi+phi2(i)-pi;
end
for i=1:num2
    phi2(i)=pi*3/4+pi*3/4-phi2(i);
end
F_D=PofD*A_piston;
F_Ds=F_D*exp(1i*(-1)*pi/4);
for i=1:num2
    rA = OA*exp(1i*phi2(i));
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
    Mcs(i)=cdot(F_Ds(i),vD)/omg1;
end
deg2_8pi = zeros(num2*4,1);
Mc_8pi = zeros(num2*4,1);
for i=1:num2
    deg2_8pi(i) = phi2(i);
    deg2_8pi(i+num2) = 2*pi+phi2(i);
    deg2_8pi(i+num2*2) = 4*pi+phi2(i);
    deg2_8pi(i+num2*3) = 6*pi+phi2(i);
    Mc_8pi(i) = Mcs(i);
    Mc_8pi(i+num2) = Mcs(i);
    Mc_8pi(i+num2*2) = Mcs(i);
    Mc_8pi(i+num2*3) = Mcs(i);
end
[deg2_fix, index] = unique(deg2_8pi-2*pi); 
Mc_unique=interp1(deg2_fix,Mc_8pi(index),deg1_unique','pchip');
Ac=cumtrapz(deg1_unique(2),Mc_unique);
% % plot(deg1_unique,Ac);
% %h1 = length of last point in Ac diagram
% %h2 = length of last point in Ad diagram
% % TLX_P_D_real= h2/h1*TLX_P_D_fake;
h1=Ac(length(Ac));
h2=Ad(length(Ad));
PD_real=h2*PofD/h1;
FD_real=PD_real*A_piston;
FD_real_8pi=[FD_real;FD_real;FD_real;FD_real];
FD_unique=interp1(deg2_fix,FD_real_8pi(index),deg1_unique','pchip');
Mc_real=Mc_unique*(h2/h1);
Ac_real=cumtrapz(deg1_unique(2),Mc_real);
Ac_deg1 = zeros(length(deg1_unique),2);
for i=1:length(deg1_unique)
    Ac_deg1(i,1) = deg1_unique(i);
    Ac_deg1(i,2) = Ac_real(i);
end
FD_deg1 = zeros (length(deg1_unique),2);
for i=1:length(deg1_unique)
    FD_deg1(i,1)=deg1_unique(i);
    FD_deg1(i,2)=FD_unique(i);
end

% %For energy calculation
Js_unique= zeros (length(deg1_unique),1);
[deg1_Js, index] = unique(Js4pi(:,1)*(4*pi)/720); 
for i=1:length(deg1_unique)
    Js_unique(i) = interp1(deg1_Js,Js4pi(index,2),deg1_unique(i),'pchip');
end
E=Ad_deg1(:,2)-Ac_deg1(:,2) +1/2*omg1^2*Js4pi(1,2);
% %chon delta(he so khong deu cho phep)=1/80
% % Chon muy J =1; muy E=1;
psi_max=(atan((muyJ/(2*muyE)*pi/180 ))*omg1^2*(1+delta/2)^2)*pi/180;
psi_min=(atan((muyJ/(2*muyE)*pi/180 ))*omg1^2*(1-delta/2)^2)*pi/180;
% %Find the right tangent
% %y=ax+b with a=tan(psi) and b is the point on E axis 
% % (a is psi_maxor and b is psi_min)
pointA=E-tan(psi_max)*Js_unique;
pointB=E-tan(psi_min)*Js_unique;
a=max(pointA);
b=min(pointB);
x_td = linspace(-0.1,0.15,26);
ytmax=tan(psi_max)*x_td+a;
ytmin=tan(psi_min)*x_td+b;
Jd=muyJ*(a-b)/(tan(psi_max)-tan(psi_min))-Js4pi(1,2);

% % FOR PB_s and PD_s
% % plot(s_B,P_B);
% % plot(s_D,PofD);
% 
% % FOR FB and FD
% % plot(deg1_unique,FB_unique);
% % plot(deg1_unique,FD_unique);
% 
% % FOR MOMENT DONG, MOMENT CAN
% % plot(deg1_unique,Md_unique);
% % plot(deg1_unique,Mc_unique);
% 
% % FOR CONG DONG, CONG CAN
% % plot(Ad_deg1(:,1),Ad_deg1(:,2),'b',Ac_deg1(:,1),Ac_deg1(:,2),'r');
% 
% % FOR E-J and flywheel
% % plot(deg1_unique,E,'k');
% grid on, hold on
% plot(Js_unique,E,'k');
% plot(x_td,ytmax,'r',x_td,ytmin,'r');

function P = cdot(F, v)
  P = real(F)*real(v) + imag(v)*imag(F);
end