function formulate(omg1, OA, AC, AB, BC)

    syms x y phi(t) t CD theta
    
    params = readtable('parameters.txt');
    
    alpha = deg2rad(params{1,1});
%     g = params{1,11};
%     m1 = params{1,6};         J1 = params{1,12};    Q1 = -m1*g*[0,1,0];
%     m2 = params{1,7};         J2 = params{1,13};    Q2 = -m2*g*[0,1,0];
%     m3 = params{1,8};         J3 = params{1,14};    Q3 = -m3*g*[0,1,0];
%     m4 = params{1,9};         J4 = params{1,15};    Q4 = -m4*g*[0,1,0];
%     m5 = params{1,10};        J5 = params{1,16};    Q5 = -m5*g*[0,1,0];
    
    xOB = pi/2 - alpha/2;
    xOD = pi/2 + alpha/2;
    
    old = {diff(phi, t, 2), diff(phi), phi};
    new = {0, omg1, theta};
    
    rA = OA*[cos(phi(t)), sin(phi(t)), 0];
    vA = diff(rA, t);
    aA = diff(vA, t);
    fpA = matlabFunction(subs(rA, old, new));
    fvA = matlabFunction(subs(vA, old, new));
    faA = matlabFunction(subs(aA, old, new));
    
    rB = [x, y, 0];
    eq1 = (rA(1)-rB(1))^2 + (rA(2)-rB(2))^2 == AB^2;
    eq2 = y/x == tan(xOB);
    sol = solve([eq1, eq2], [x, y]);
    rB = simplify(subs(rB, [x, y], [sol.x(2), sol.y(2)]));
    vB = diff(rB, t);
    aB = diff(vB, t);
    fpB = matlabFunction(subs(rB, old, new));
    fvB = matlabFunction(subs(vB, old, new));
    faB = matlabFunction(subs(aB, old, new));
    
    rC = [x, y, 0];
    eq1 = (rA(1)-rC(1))^2 + (rA(2)-rC(2))^2 == AC^2;
    eq2 = (rB(1)-rC(1))^2 + (rB(2)-rC(2))^2 == BC^2;
    sol = solve([eq1, eq2], [x, y]);
    rC = simplify(subs(rC, [x, y], [sol.x(2), sol.y(2)]));
    vC = diff(rC, t);
    aC = diff(vC, t);
    fpC = matlabFunction(subs(rC, old, new));
    fvC = matlabFunction(subs(vC, old, new));
    faC = matlabFunction(subs(aC, old, new));
    
    rD = [x, y, 0];
    eq1 = (rC(1)-rD(1))^2 + (rC(2)-rD(2))^2 == CD^2;
    eq2 = y/x == tan(xOD);
    sol = solve([eq1, eq2], [x, y]);
    rD = simplify(subs(rD, [x, y], [sol.x(1), sol.y(1)]));
    vD = diff(rD, t);
    aD = diff(vD, t);
    fpD = matlabFunction(subs(rD, old, new), 'Vars', [theta, CD]);
    fvD = matlabFunction(subs(vD, old, new), 'Vars', [theta, CD]);
    faD = matlabFunction(subs(aD, old, new), 'Vars', [theta, CD]);
    
    omg2 = [0, 0, x];
    eq = vB == cross(omg2, rB-rA) + vA;
    sol = solve(eq, x);
    omg2 = subs(omg2, x, sol);
    alp2 = diff(omg2, t);
    fomg2 = matlabFunction(subs(omg2, old, new));
    falp2 = matlabFunction(subs(alp2, old, new));
    
    omg4 = [0, 0, x];
    eq = vD == cross(omg4, rD-rC) + vC;
    sol = solve(eq, x);
    omg4 = subs(omg4, x, sol);
    alp4 = diff(omg4, t);
    fomg4 = matlabFunction(subs(omg4, old, new), 'Vars', [theta, CD]);
    falp4 = matlabFunction(subs(alp4, old, new), 'Vars', [theta, CD]);

    save('./outputs/functions.mat', 'fpA', 'fpB', 'fpC', 'fpD',...
        'fvA', 'fvB', 'fvC', 'fvD', 'faA', 'faB', 'faC', 'faD',...
        'fomg2', 'fomg4', 'falp2', 'falp4');
end