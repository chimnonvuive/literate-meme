function [aB, aC, aE, alp2, alp3] = d_acc(phi)
[vB, vC, vE, omg2, omg3] = d_vel(phi);

% find aB, aC, aE
aB = diff(vB); aC = diff(vC); aE = diff(vE);

% find alp2, alp3
alp2 = diff(omg2); alp3 = diff(omg3);
end
