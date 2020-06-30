global AB BC CD CE rD omg1 alp1
AB = .15; BC = .35; CD = .3; CE = .15;
rD = [.3,.3,0]; omg1 = [0,0,1]; alp1 = [0,0,0];
% 
% i=0;theta1s = zeros(201,1), aCs = theta1s, vCs = aCs;
% 
% for theta1=0:(2*pi/200):2*pi
%     i=i+1;
%     [~,vC,~,~,~] = vel(theta1); vCs(i,1) = vC(1);
%     [~,aC,~,~,~] = acc(theta1); aCs(i,1) = aC(1);
%     theta1s(i)=theta1;
% end
% plot(theta1s, vCs(:,1))
% plot(theta1s, aCs(:,1))
% 
% for theta1=0:2*pi/32:2*pi
%     hold on, axis equal
%     [rB, rC, rE] = pos(theta1);
%     plot([0, rB(1)], [0, rB(2)], 'r-o');
%     plot([rB(1), rC(1)], [rB(2), rC(2)], 'g-o');
%     plot([rD(1), rE(1)], [rD(2), rE(2)], 'b-o');
% end
% 
% hold on, axis equal
% [rB, rC, rE] = pos(pi/4);
% plot([0, rB(1)], [0,rB(2)], 'r-o')
% plot([rB(1), rC(1)], [rB(2), rC(2)], 'g-o')
% plot([rD(1), rE(1)], [rD(2), rE(2)], 'b-o')
% 
syms theta(t) t
theta1 = pi/4;
old = {diff(theta, t, 2), diff(theta), theta}; new = {alp1(3), omg1(3), theta1};
% 
% [rB, rC, rE] = pos(theta1)
% [d_rB, d_rC, d_rE] = d_pos(theta1);
% d_rB = double(subs(d_rB, theta, theta1))
% d_rC = double(subs(d_rC, theta, theta1))
% d_rE = double(subs(d_rE, theta, theta1))
% 
% [vB, vC, vE, omg2, omg3] = vel(theta1)
% [c_vB, c_vC, c_vE, c_omg2, c_omg3] = c_vel(theta1)
% [d_vB, d_vC, d_vE, d_omg2, d_omg3] = d_vel(theta1);
% d_vB = double(subs(d_vB, old, new))
% d_vC = double(subs(d_vC, old, new))
% d_vE = double(subs(d_vE, old, new))
% d_omg2 = double(subs(d_omg2, old, new))
% d_omg3 = double(subs(d_omg3, old, new))
% 
[aB, aC, aE, alp2, alp3] = acc(theta1)
[c_aB, c_aC, c_aE, c_alp2, c_alp3] = c_acc(theta1)
[d_aB, d_aC, d_aE, d_alp2, d_alp3] = d_acc(theta1);
d_aB = double(subs(d_aB, old, new))
d_aC = double(subs(d_aC, old, new))
d_aE = double(subs(d_aE, old, new))
d_alp2 = double(subs(d_alp2, old, new))
d_alp3 = double(subs(d_alp3, old, new))


function [rB, rC, rE] = pos(theta1)
global AB BC CD CE rD, syms x y
% set variables
rB = AB*[cos(theta1), sin(theta1), 0];
rC = [x, y, 0]; rE = [x, y, 0];
% find rC
eq1 = (x-rB(1))^2 + (y-rB(2))^2 == BC^2;
eq2 = (x-rD(1))^2 + (y-rD(2))^2 == CD^2;
sol = solve([eq1,eq2], [x,y]);
% return the values
if 0>sol.x(1)
    rC = double(subs(rC, [x, y], [sol.x(1), sol.y(1)]));
else
    rC = double(subs(rC, [x, y], [sol.x(2), sol.y(2)]));    
end
% find rE
eq1 = (rD(2)-rE(2))/(rD(1)-rE(1)) == (rD(2)-rC(2))/(rD(1)-rC(1));
eq2 = (rC(1)-rE(1))^2 + (rC(2)-rE(2))^2 == CE^2;
sol = solve([eq1,eq2], [x,y]);
% return the values
if rC(1)>double(sol.x(1))
    rE = double(subs(rE, [x, y], [sol.x(1), sol.y(1)]));
else
    rE = double(subs(rE, [x, y], [sol.x(2), sol.y(2)]));
end
end
function [vB, vC, vE, omg2, omg3] = vel(theta1)
global omg1 rD, syms x y
% set variables
[rB, rC, rE] = pos(theta1);
omg2 = [0,0,x]; omg3 = [0,0,y];
vB = cross(omg1, rB);
vC = cross(omg3, rC-rD);
vE = cross(omg3, rE-rD);
% find omg2, omg3
eq = vB + cross(omg2, rC-rB) == vC;
sol = solve([eq(1),eq(2)], [x,y]);
% return the values
omg2 = double(subs(omg2, sol.x));
omg3 = double(subs(omg3, sol.y));
vC = double(subs(vC, sol.y));
vE = double(subs(vE, sol.y));
end
function [aB, aC, aE, alp2, alp3] = acc(theta1)
global omg1 rD alp1, syms x y
% set variables
[rB, rC, rE] = pos(theta1);
[~,~,~, omg2, omg3] = vel(theta1);
alp2 = [0,0,x]; alp3 = [0,0,y];
aB = cross(alp1, rB) - omg1(3)^2*rB;
aC = cross(alp3, rC-rD) - omg3(3)^2*(rC-rD);
aE = cross(alp3, rE-rD) - omg3(3)^2*(rE-rD);
% find alp2, alp3
eq = aB + cross(alp2, rC-rB) - omg2(3)^2*(rC-rB) == aC;
sol = solve([eq(1),eq(2)], [x,y]);
% find aC, aE
aC = double(subs(aC, sol.y));
aE = double(subs(aE, sol.y));
alp2 = double(subs(alp2, sol.x));
alp3 = double(subs(alp3, sol.y));
end


function [rB, rC, rE] = d_pos(theta1)
global AB BC CD CE rD, syms x y theta(t)
% set variables
rB = AB*[cos(theta(t)), sin(theta(t)), 0];
rC = [x, y, 0]; rE = [x, y, 0];
% find rC
eq1 = (rB(1)-x)^2 + (rB(2)-y)^2 == BC^2;
eq2 = (x-rD(1))^2 + (y-rD(2))^2 == CD^2;
sol = solve([eq1,eq2], [x,y]);
if 0>subs(sol.x(1), theta, theta1)
    rC = subs(rC, [x,y], [sol.x(1),sol.y(1)]);
else
    rC = subs(rC, [x,y], [sol.x(2),sol.y(2)]);
end
% find rE
eq1 = (rD(2)-y)/(rD(1)-x) == (rD(2)-rC(2))/(rD(1)-rC(1));
eq2 = (rC(1)-x)^2 + (rC(2)-y)^2 == CE^2;
sol = solve([eq1,eq2], [x,y]);
if subs(rC(1), theta, theta1) > subs(sol.x(1), theta, theta1)
    rE = subs(rE, [x, y], [sol.x(1), sol.y(1)]);
else
    rE = subs(rE, [x, y], [sol.x(2), sol.y(2)]);
end
end
function [vB, vC, vE, omg2, omg3] = d_vel(theta1)
global rD
[rB, rC, rE] = d_pos(theta1);
theta2 = atan((rC(2)-rB(2))/(rC(1)-rB(1)));
theta3 = atan((rC(2)-rD(2))/(rC(1)-rD(1)));
% find vB, vC, vE, omg2, omg3
vB = diff(rB); vC = diff(rC); vE = diff(rE);
omg2 = [0,0,diff(theta2)]; omg3 = [0,0,diff(theta3)];
end
function [aB, aC, aE, alp2, alp3] = d_acc(theta1)
[vB, vC, vE, omg2, omg3] = d_vel(theta1);
% find aB, aC, aE, alp2, alp3
aB = diff(vB); aC = diff(vC); aE = diff(vE);
alp2 = diff(omg2); alp3 = diff(omg3);
end


function [vB, vC, vE, omg2, omg3] = c_vel(theta1)
global rD omg1, syms x y z
% set variables
[rB, rC, rE] = pos(theta1); vB = cross(omg1, rB);
omg12 = [0,0,x]; omg23 = [0,0,y]; omg30 = [0,0,z];
% find omg12, omg23, omg30
eq1 = omg1 + omg12 + omg23 + omg30 == 0;
eq2 = cross(rB, omg12) + cross(rC, omg23) + cross(rD, omg30) == 0;
sol = solve([eq1(3),eq2(1),eq2(2)], [x,y,z]);
% return the values
omg12 = double(subs(omg12, sol.x));
omg23 = double(subs(omg23, sol.y));
omg30 = double(subs(omg30, sol.z));
omg2 = omg1 + omg12; omg3 = -omg30;
vC = cross(omg3, rC-rD);
vE = cross(omg3, rE-rD);
end
function [aB, aC, aE, alp2, alp3] = c_acc(theta1)
global omg1 alp1 rD, syms x y z
% set variables
[rB, rC, rE] = pos(theta1);
[~,~,~, omg2, omg3] = vel(theta1);
aB = cross(alp1, rB) - omg1(3)^2*rB;
alp12 = [0,0,x]; alp23 = [0,0,y]; alp30 = [0,0,z];
% find alp12, alp23, alp30
eq1 = alp1 + alp12 + alp23 + alp30 == 0;
eq2 = cross(rB, alp12) + cross(rC, alp23) + cross(rD, alp30) - omg1(3)^2*rB - omg2(3)^2*(rC-rB) - omg3(3)^2*(rD-rC) == 0;
sol = solve([eq1(3),eq2(1),eq2(2)], [x,y,z]);
% return the values
alp12 = double(subs(alp12, sol.x));
alp23 = double(subs(alp23, sol.y));
alp30 = double(subs(alp30, sol.z));
alp2 = alp1 + alp12; alp3 = -alp30;
aC = cross(alp3, rC-rD) - omg3(3)^2*(rC-rD);
aE = cross(alp3, rE-rD) - omg3(3)^2*(rE-rD);
end
