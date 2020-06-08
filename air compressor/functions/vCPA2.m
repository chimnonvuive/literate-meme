function [v1_T, omg2] = vCPA2(v_known, r1, r2, v2_T, omg1)

%     Complex-algebraic velocity analysis, Theory of Machines and Mechanisms
%     The unknowns are v1_T, omg2
%     v_known + v_r1 = v_r2
    
    A = [-cos(angle(r1)), -imag(r2);
         -sin(angle(r1)),  real(r2)];
    b = [real(v_known) - ( imag(r1)*omg1 + cos(angle(r2))*v2_T);
         imag(v_known) - (-real(r1)*omg1 + sin(angle(r2))*v2_T)];
    res = linsolve(A, b); v1_T = res(1); omg2 = res(2);
end