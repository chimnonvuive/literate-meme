%Energy calculation
HB = 39; alpha = 90*pi/180; beta = 90*pi/180; n1 = 500;
iter = 20000;
AB = 104.75; 
CD = 85.653; HD = 39.0362; % iter = 20000
OA = HB/2; AC = OA;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;

load('Ac_energy.mat', 'Ac_energy');
load('Ad_energy.mat', 'Ad_energy');
load('Js_0_to_720.mat', 'Js_0_to_720');
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
%Chon muy J =10; if draw on a square, then y(real)=1 and x(real)=1
delta=1/80;
muyE=1;
muyJ=1;
psi_max=atan(muyJ/(2*muyE))*omg1^2*(1+delta/2)^2;
psi_min=atan(muyJ/(2*muyE))*omg1^2*(1-delta/2)^2;

%Find the right tangent
x_td = zeros (26,1);
for i=1:26
    x_td(i)=-0.1+(i-1)*0.01;
end
[Emax,indmax]=max(E);
[Emin,indmin]=min(E);
ytmax=tan(psi_max)*(x_td-new_Js(indmax))+Emax;
ytmin=tan(psi_min)*(x_td-new_Js(indmin))+Emin;
plot(x_td,ytmax,'r',x_td,ytmin,'r');
grid on
Jd=muyJ*(ytmax(11)-ytmin(11))/(tan(psi_max)-tan(psi_min))-Js_0_to_720(1,2);