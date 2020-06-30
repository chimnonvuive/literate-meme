function [B_angle_1, A_1, B_angle_2, A_2] = CPA2(A_angle, B, c)
  C = abs(c); C_angle = angle(c);
  asin_1 =     asin(C*sin(C_angle-A_angle)/B);
  asin_2 =  pi-asin(C*sin(C_angle-A_angle)/B);
  B_angle_1 = A_angle + asin_1;
  B_angle_2 = A_angle + asin_2;
  A_1 = C*cos(C_angle-A_angle) - B*cos(B_angle_1-A_angle);
  A_2 = C*cos(C_angle-A_angle) - B*cos(B_angle_2-A_angle);
end