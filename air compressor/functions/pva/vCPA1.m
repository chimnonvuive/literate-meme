function [v1_T, omg1] = vCPA1(v_known, r1, r2, v2_T, omg2)
    
%     Complex-algebraic velocity analysis, Theory of Machines and Mechanisms
%     The unknowns are v1_T, omg1
%     v_known + v_r1 = v_r2

    A = [-cos(angle(r1)),  imag(r1);
         -sin(angle(r1)), -real(r1)];
    b = [real(v_known) - (cos(angle(r2))*v2_T - imag(r2)*omg2); 
         imag(v_known) - (sin(angle(r2))*v2_T + real(r2)*omg2)];
    res = linsolve(A, b); v1_T = res(1); omg1 = res(2);
end