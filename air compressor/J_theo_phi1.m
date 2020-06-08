%J theo phi1
HB = 39; alpha = 90*pi/180; beta = 90*pi/180; n1 = 500; iter = 20000;
AB = 104.75; 
CD = 85.653; HD = 39.0362; % iter = 20000
OA = HB/2; AC = OA;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;
m2=0.5; %kg, g=10 m/s^2
m3=m2; m4=m2; m5=m2;
J1 = 0.12; %Kg.m^2
J2 = J1; J3 = J1; J4 = J1; J5 = J1;
AB=AB/1000; HB=HB/1000; CD=CD/1000; HD=HD/1000;OA=OA/1000;AC=AC/1000;
BC=BC/1000;
Js = zeros(iter,1);

i=1;
for theta=linspace(0,2*pi,iter+1)
    rA = OA*exp(1i*theta);
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
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
hold on
plot(linspace(0,360,iter+1),Js,'k');
plot(linspace(360,720,iter+1),Js,'k');