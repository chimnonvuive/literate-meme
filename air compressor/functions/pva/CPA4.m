function [A_angle_1, B_angle_1, A_angle_2, B_angle_2] = CPA4(A, B, c)
    
%     Complex Polar Algebra, method 4, Theory of Machines and Mechanisms
%     The unknowns are A_angle, B_angle
%     Using this case leads to 2 different results. Be sure to check out for
%     appropriate conditions to come up with the right solution.
    
    C = abs(c); C_angle = angle(c); pm = [1, -1];
    B_angle = C_angle + pm*acos((C^2+B^2-A^2)/(2*B*C));
    A_angle = C_angle + pm*acos((C^2+A^2-B^2)/(2*A*C));
    if A*sin(A_angle(1) - C_angle) == -B*sin(B_angle(1) - C_angle)
        A_angle_1 = A_angle(1); B_angle_1 = B_angle(1);
        A_angle_2 = A_angle(2); B_angle_2 = B_angle(2);
    else
        A_angle_1 = A_angle(1); B_angle_1 = B_angle(2);
        A_angle_2 = A_angle(2); B_angle_2 = B_angle(1);
    end
end