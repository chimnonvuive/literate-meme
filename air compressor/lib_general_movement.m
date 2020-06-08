addpath('./functions')

HB = 39; alpha = 90*pi/180; beta = 90*pi/180; n1 = 500; iter = 100;
AB = 104.395; CD = 94.8311; HD = 39.0074;
OA = HB/2; AC = OA/2;  BC = sqrt(AC^2+AB^2-2*AC*AB*cos(beta));
xOB = pi/2 - alpha/2; xOD = pi/2 + alpha/2; omg1 = n1*pi/30;

figure, axis equal
sys(iter) = struct('cdata',[],'colormap',[]);

v = VideoWriter('./outputs/movements.avi'); v.Quality = 100; open(v)
phis = linspace(0,360,100);

for i = 1:length(phis)
    rA = OA*exp(1i*phis(i)*pi/180);
    [~,~,~,OB] = CPA2(xOB, AB, rA);
    rB = OB*exp(1i*xOB);
    [xAC,~,~,~] = CPA4(AC, BC, rB-rA);
    rC = rA + AC*exp(1i*xAC);
    [~,~,~,OD] = CPA2(xOD, CD, rC);
    rD = OD*exp(1i*xOD);
    
%     ----- plot position -----
    
    link1 = [0,rA]; 
    link2_AB = [rA,rB];
    link2_AC = [rA,rC];
    link2_BC = [rB,rC]; 
    link4 = [rC,rD];
%     
    hold on, axis equal
    plot(real(link1), imag(link1), 'bo-')
    plot(real(link2_AB), imag(link2_AB), 'ro-')
    plot(real(link2_AC), imag(link2_AC), 'ro-')
    plot(real(link2_BC), imag(link2_BC), 'ro-')
    plot(real(link4), imag(link4), 'ko-')
    xlim([-150,150]); ylim([-150,150]);
    drawnow; sys(i) = getframe; clf
end
writeVideo(v,sys); close(v)