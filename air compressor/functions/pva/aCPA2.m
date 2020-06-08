function [a1_T, alp2] = aCPA2(a_known, r1, r2, v1_T, omg1, v2_T, omg2, a2_T, alp1)
  a1_known = centripetal(r1, omg1) + coriolis(r1, v1_T, omg1);
  a2_known = centripetal(r2, omg2) + coriolis(r2, v2_T, omg2);
  free_coeff = a_known + a1_known - a2_known;
  
  A = [-cos(angle(r1)), -imag(r2);
       -sin(angle(r1)),  real(r2)];
  b = [real(free_coeff) - ( imag(r1)*alp1 + cos(angle(r2))*a2_T);
       imag(free_coeff) - (-real(r1)*alp1 + sin(angle(r2))*a2_T)];
  res = linsolve(A, b); a1_T = res(1); alp2 = res(2);
end
