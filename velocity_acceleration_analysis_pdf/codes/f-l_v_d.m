function [vB, vC, vE, omg2, omg3] = d_vel(phi)
global rD

[rB, rC, rE] = d_pos(phi);

% find vB, vC, vE
vB = diff(rB); vC = diff(rC); vE = diff(rE);

% find omg2, omg3
phi2 = atan((rC(2)-rB(2))/(rC(1)-rB(1)));
phi3 = atan((rC(2)-rD(2))/(rC(1)-rD(1)));
omg2 = [0,0,diff(phi2)];
omg3 = [0,0,diff(phi3)];
end
