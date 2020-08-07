function formulate(omg1, alpha, OA, AC, AB, BC, m2)

    syms x y phi(t) t CD theta
    
    xOB = pi/2 - alpha/2;
    xOD = pi/2 + alpha/2;
    
    mAC = AC/(AB+AC)*m2;
    mAB = AB/(AB+AC)*m2;
    
    old = {diff(phi, t, 2), diff(phi), phi};
    new = {0, omg1, theta};
    
    rA = OA*[cos(phi(t)), sin(phi(t)), 0];
    vA = diff(rA);
    aA = diff(vA);
    fpA = matlabFunction(subs(rA, old, new));
    fvA = matlabFunction(subs(vA, old, new));
    faA = matlabFunction(subs(aA, old, new));
    
    rB = [x, y, 0];
    eq1 = (rA(1)-rB(1))^2 + (rA(2)-rB(2))^2 == AB^2;
    eq2 = y/x == tan(xOB);
    sol = solve([eq1, eq2], [x, y]);
    rB = simplify(subs(rB, [x, y], [sol.x(2), sol.y(2)]));
    vB = diff(rB);
    aB = diff(vB);
    fpB = matlabFunction(subs(rB, old, new));
    fvB = matlabFunction(subs(vB, old, new));
    faB = matlabFunction(subs(aB, old, new));
    
    rC = [x, y, 0];
    eq1 = (rA(1)-rC(1))^2 + (rA(2)-rC(2))^2 == AC^2;
    eq2 = (rB(1)-rC(1))^2 + (rB(2)-rC(2))^2 == BC^2;
    sol = solve([eq1, eq2], [x, y]);
    rC = simplify(subs(rC, [x, y], [sol.x(2), sol.y(2)]));
    vC = diff(rC);
    aC = diff(vC);
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
    
    rAC = simplify(rA-rC); rAB = simplify(rA-rB);
    p2 = (mAC*rAC + mAB*rAB)/m2;
    phi2 = simplify(atan((p2(2)-rA(2))/(p2(1)-rA(1))));
    omg2 = diff(phi2);
    alp2 = diff(omg2);
    fomg2 = matlabFunction(subs(omg2, old, new));
    falp2 = matlabFunction(subs(alp2, old, new));
    p4 = (rC+rD)/2;
    phi4 = simplify(atan((p4(2)-rC(2))/(p4(1)-rC(1))));
    omg4 = diff(phi4, t);
    alp4 = diff(omg4, t);
    fomg4 = matlabFunction(subs(omg4, old, new), 'Vars', [theta, CD]);
    falp4 = matlabFunction(subs(alp4, old, new), 'Vars', [theta, CD]);
    
    save('./outputs/functions.mat', 'fpA', 'fpB', 'fpC', 'fpD',...
        'fvA', 'fvB', 'fvC', 'fvD', 'faA', 'faB', 'faC', 'faD',...
        'fomg2', 'fomg4', 'falp2', 'falp4');
end