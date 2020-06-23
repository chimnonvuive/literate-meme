global AB BC omg1 alp1
AB = .5; BC = 1; omg1 = [0,0,1]; alp1 = [0,0,-1];

% hold on, axis equal
% [rB, rC] = pos(pi/3);
% plot([0, rB(1)], [0,rB(2)], 'r-o')
% plot([rB(1), rC(1)], [rB(2), rC(2)], 'b-o');

% for phi1=0:2*pi/32:2*pi
%     hold on, axis equal
%     [rB, rC] = pos(phi1);
%     plot([0, rB(1)], [0,rB(2)], 'r-o');
%     plot([rB(1), rC(1)], [rB(2), rC(2)], 'b-o');
% end

% i=0; phi1s=zeros(201); vCs=phi1s; aCs=phi1s;
% for phi1=0:(2*pi/200):2*pi
%     i=i+1;
%     [~,vC,~] = vel(phi1); vCs(i) = vC(1);
%     [~,aC,~] = acc(phi1); aCs(i) = aC(1);
%     phi1s(i)=phi1;
% end
% plot(phi1s, vCs)
% plot(phi1s, aCs)
% 
% syms phi(t)
% phi1 = pi/3;
% old = {diff(phi, t, 2), diff(phi), phi}; new = {alp1(3), omg1(3), phi1};
% 
% [rB, rC] = pos(phi1)
% [d_rB, d_rC] = d_pos(phi1);
% d_rB = double(subs(d_rB, old, new))
% d_rC = double(subs(d_rC, old, new))
% 
% [vB, vC, omg2] = vel(phi1)
% [c_vB, c_vC, c_omg2] = c_vel(phi1)
% [d_vB, d_vC, d_omg2] = d_vel(phi1);
% d_vB = double(subs(d_vB, old, new))
% d_vC = double(subs(d_vC, old, new))
% d_omg2 = double(subs(d_omg2, old, new))
% 
% [aB, aC, alp2] = acc(phi1)
% [c_aB, c_aC, c_alp2] = c_acc(phi1)
% [d_aB, d_aC, d_alp2] = d_acc(phi1);
% d_aB = double(subs(d_aB, old, new))
% d_aC = double(subs(d_aC, old, new))
% d_alp2 = double(subs(d_alp2, old, new))

function [rB, rC] = pos(phi1)
global AB BC, syms x
% set variables
rB = AB*[cos(phi1), sin(phi1), 0]; rC = [x, 0, 0];
% find rC
eq = (rB(1)-rC(1))^2 + rB(2)^2 == BC^2;
sol = solve(eq, x);
% return the values
if subs(rC(1), sol(1)) > subs(rC(1), sol(2))
    rC = double(subs(rC, sol(1)));
else
    rC = double(subs(rC, sol(2)));
end
end
function [vB, vC, omg2] = vel(phi1)
global omg1, syms x y
% set variables
[rB, rC] = pos(phi1);
vB = cross(omg1, rB);
vC = [x,0,0]; omg2 = [0,0,y];
% find vC, omg2
eq = vB + cross(omg2, rC-rB) == vC;
sol = solve([eq(1),eq(2)], [x,y]);
% return the values
vC = double(subs(vC, sol.x));
omg2 = double(subs(omg2, sol.y));
end
function [aB, aC, alp2] = acc(phi1)
global omg1 alp1, syms x y
% set variables
[rB, rC] = pos(phi1);
[~,~, omg2] = vel(phi1);
aB = cross(alp1, rB) - norm(omg1)^2*rB;
aC = [x,0,0]; alp2 = [0,0,y];
% find aC, alp2
eq = aB + cross(alp2, rC-rB) - omg2(3)^2*(rC-rB) == aC;
sol = solve(eq(1), eq(2), x, y);
% return the values
aC = double(subs(aC, sol.x));
alp2 = double(subs(alp2, sol.y));
end


function [rB, rC] = d_pos(phi1)
global AB BC, syms x phi(t)
% set variables
rB = AB*[cos(phi(t)), sin(phi(t)), 0];
rC = [x,0,0];
% find rC
eq = (rB(1)-x)^2 + rB(2)^2 == BC^2;
sol = solve(eq, x);
if subs(sol(1),phi,phi1) > subs(sol(2),phi,phi1)
    rC = subs(rC, sol(1));
else
    rC = subs(rC, sol(2));
end
end
function [vB, vC, omg2] = d_vel(phi1)
[rB, rC] = d_pos(phi1); phi2 = atan((rC(2)-rB(2))/(rC(1)-rB(1)));
% find vB, vC, omg2
vB = diff(rB); vC = diff(rC); omg2 = [0,0,diff(phi2)];
end
function [aB, aC, alp2] = d_acc(phi1)
[vB, vC, omg2] = d_vel(phi1);
% find aB, aC, alp2
aB = diff(vB); aC = diff(vC); alp2 = diff(omg2);
end


function [vB, vC, omg2, v23] = c_vel(phi1)
global omg1, syms x y z
% set variables
[rB, rC] = pos(phi1); vB = cross(omg1, rB);
omg12 = [0,0,x]; omg23 = [0,0,y]; v23 = [z,0,0];
% find omg12, omg23, v23
eq1 = omg1 + omg12 + omg23 == 0;
eq2 = cross(rB, omg12) + cross(rC, omg23) + v23 == 0;
sol = solve([eq1(3),eq2(1),eq2(2)], [x,y,z]);
% return the values
omg12 = double(subs(omg12, sol.x));
omg23 = double(subs(omg23, sol.y));
v23 = double(subs(v23, sol.z));
vC = -v23;
omg2 = omg1 + omg12;
end
function [aB, aC, alp2] = c_acc(phi1)
global omg1 alp1, syms x y z
% set variables
[rB, rC] = pos(phi1);
[~,~,omg2,~] = c_vel(phi1);
aB = cross(alp1,rB) - omg1(3)^2*rB;
alp12 = [0,0,x]; alp23 = [0,0,y]; a23 = [z,0,0];
% find alp12, alp23, a23
eq1 = alp1 + alp12 + alp23 == 0;
eq2 = cross(rB,alp12) + cross(rC,alp23) + a23 - omg1(3)^2*rB - omg2(3)^2*(rC-rB) == 0; 
sol = solve([eq1(3),eq2(1),eq2(2)], [x,y,z]);
% return the values
alp12 = double(subs(alp12, sol.x));
alp23 = double(subs(alp23, sol.y));
a23 = double(subs(a23, sol.z));
aC = -a23;
alp2 = alp1 + alp12;
end
