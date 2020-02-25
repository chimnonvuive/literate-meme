function [vB, vC, omg2] = d_vel(phi)
[rB, rC] = d_pos(phi);
phi2 = atan((rC(2)-rB(2))/(rC(1)-rB(1)));
% find vB, vC, omg2
vB = diff(rB); vC = diff(rC);
omg2 = [0,0,diff(phi2)]; 
end


function [aB, aC, alp2] = d_acc(phi)
[vB, vC, omg2] = d_vel(phi);
% find aB, aC, alp2
aB = diff(vB); aC = diff(vC); alp2 = diff(omg2);
end
