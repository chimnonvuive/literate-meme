% [fA,fB,fC,fD] = dick(pi/2, 0.0195, 0.0195, 0.10475, 0.096489312879717, 0.085653);
% rA = zeros(100,3);
% rB = zeros(100,3);
% rC = zeros(100,3);
% rD = zeros(100,3);
% phi = linspace(0,2*pi,100);
% iA = zeros(100,1);
% iB = zeros(100,1);
% iC = zeros(100,1);
% iD = zeros(100,1);
% 
% for i=1:100
%     rA(i,:) = fA(phi(i)); iA(i) = rA(i,1) + 1j*rA(i,2);
%     rB(i,:) = fB(phi(i)); iB(i) = rB(i,1) + 1j*rB(i,2);
%     rC(i,:) = fC(phi(i)); iC(i) = rC(i,1) + 1j*rC(i,2);
%     rD(i,:) = fD(phi(i)); iD(i) = rD(i,1) + 1j*rD(i,2);
% end
% 
% animation(iA,iB,iC,iD);
% 
% 
function formulate(alpha, OA, AC, AB, BC)
    syms x y phi CD
    xOB = pi/2 - alpha/2;
    xOD = pi/2 + alpha/2;
    
    rA = OA*[cos(phi), sin(phi), 0];
    fpA = matlabFunction(rA);
    fvA = matlabFunction(diff(rA));
    faA = matlabFunction(diff(rA, phi, 2));
    
    rB = [x, y, 0];
    eq1 = (rA(1)-rB(1))^2 + (rA(2)-rB(2))^2 == AB^2;
    eq2 = y/x == tan(xOB);
    sol = solve([eq1, eq2], [x, y]);
    rB = subs(rB, [x, y], [sol.x(2), sol.y(2)]);
    fpB = matlabFunction(rB);
    fvB = matlabFunction(diff(rB));
    faB = matlabFunction(diff(rB, phi, 2));
    
    rC = [x, y, 0];
    eq1 = (rA(1)-rC(1))^2 + (rA(2)-rC(2))^2 == AC^2;
    eq2 = (rB(1)-rC(1))^2 + (rB(2)-rC(2))^2 == BC^2;
    sol = solve([eq1, eq2], [x, y]);
    rC = subs(rC, [x, y], [sol.x(2), sol.y(2)]);
    fpC = matlabFunction(rC);
    fvC = matlabFunction(diff(rC));
    faC = matlabFunction(diff(rC, phi, 2));
    
    rD = [x, y, 0];
    eq1 = (rC(1)-rD(1))^2 + (rC(2)-rD(2))^2 == CD^2;
    eq2 = y/x == tan(xOD);
    sol = solve([eq1, eq2], [x, y]);
    rD = subs(rD, [x, y], [sol.x(1), sol.y(1)]);
    fpD = matlabFunction(rD, 'Vars', [phi, CD]);
    fvD = matlabFunction(diff(rD), 'Vars', [phi, CD]);
    faD = matlabFunction(diff(rD, phi, 2), 'Vars', [phi, CD]);
    
    save('outputs/functions.mat', 'fpA', 'fpB', 'fpC', 'fpD',...
        'fvA', 'fvB', 'fvC', 'fvD', 'faA', 'faB', 'faC', 'faD');
end