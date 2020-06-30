global AD BC AC rC
AD = 0.35; BC = 0.20; AC = 0.15; rC = [AC,0,0];
phi1 = pi/3; omg1 = [0,0,pi]; alp1 = [0,0,0];

% hold on, axis equal
% [rB, rD] = pos(pi/3);
% plot([0, rD(1)], [0,rD(2)], 'r-o')
% plot([rB(1), rC(1)], [rB(2), rC(2)], 'b-o');

% for phi1=0:2*pi/32:2*pi
%     hold on, axis equal
%     [rB, rD] = pos(phi1);
%     plot([rC(1), rB(1)], [rC(2), rB(2)], 'b-o');
%     plot([0, rD(1)], [0, rD(2)], 'r-o');
% end

% i=0; phi1s = zeros(201, 1);
% vBs = phi1s; aBs = phi1s;
% for phi1=0:(2*pi/200):2*pi
%     i=i+1;
%     [vB2,~,~,~] = vel(phi1); vBs(i, 1) = vB2(1);
%     [aB2,~,~] = acc(phi1); aBs(i, 1) = aB2(1);
%     phi1s(i)=phi1;
% end
% plot(phi1s, vBs)
% plot(phi1s, aBs)

% syms phi(t) t
% old = {diff(phi, t, 2), diff(phi), phi}; new = {alp1(3), omg1(3), phi1};
% 
% [rB, rD] = pos(phi1)
% [d_rB, d_rD] = d_pos(phi1);
% d_rB = double(subs(d_rB, old, new))
% d_rD = double(subs(d_rD, old, new))
% 
% [vB, vD, omg3,~] = vel(phi1)
% [c_vB, c_vD, c_omg3,~] = c_vel(phi1)
% [d_vB, d_vD, d_omg3] = d_vel(phi1);
% d_vB = double(subs(d_vB, old, new))
% d_vD = double(subs(d_vD, old, new))
% omg3 = double(subs(omg3, old, new))
% 
% [aB, aD, alp3] = acc(phi1)
% [c_aB, c_aD, c_alp3] = c_acc(phi1)
% [d_aB, d_aD, d_alp3] = d_acc(phi1);
% d_aB = double(subs(d_aB, old, new))
% d_aD = double(subs(d_aD, old, new))
% d_alp3 = double(subs(d_alp3, old, new))

function [rB, rD] = pos(phi1)
global AD BC rC, syms l_AB
% set variables
rB = l_AB*[cos(phi1), sin(phi1), 0];
rD = AD*[cos(phi1), sin(phi1), 0];
% find rB
eq = (rB(1)-rC(1))^2 + (rB(2)-rC(2))^2 == BC^2;
sol = solve(eq, l_AB);
% return the values
if sol(1)>0
    rB = double(subs(rB, l_AB, sol(1)));
else
    rB = double(subs(rB, l_AB, sol(2))); end
end
function [vB, vD, omg3, vB2B1] = vel(phi1)
global rC omg1, syms x y z
% set variables
[rB, rD] = pos(phi1); vB = cross(omg3, rB-rC); vD = cross(omg1, rD);
vB2B1 = [x,y,0]; omg3 = [0,0,z];
% find vB2B1, omg3=omg2
eq1 = cross(omg3, rB-rC) == vB2B1 + cross(omg1, rB);
eq2 = dot(cross(cross(rB,[1,0,0]),rB), vB2B1) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);
% return the values
vB2B1 = double(subs(vB2B1, [x, y], [sol.x, sol.y]));
omg3 = double(subs(omg3, sol.z));
vB = double(subs(vB, sol.z));
end
function [aB, aD, alp3] = acc(phi1)
global rC alp1 omg1, syms x y z
% set variables
[rB, rD] = pos(phi1); [~,~, omg3, vB2B1] = vel(phi1);
aB2B1 = [x,y,0]; alp3 = [0,0,z];
aB1 = cross(alp1, rB) - omg1(3)^2*rB;
aD = cross(alp1, rD) - omg1(3)^2*rD;
aB2 = cross(alp3, rB-rC) - omg3(3)^2*(rB-rC);
% find aB2B1, alp3=alp2
eq1 = aB2 == aB1 + aB2B1 + 2*cross(omg1, vB2B1);
eq2 = dot(cross(cross(rB,[1,0,0]),rB), aB2B1) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);
% return the values
aB2B1 = double(subs(aB2B1, [x,y], [sol.x, sol.y]));
aB = double(subs(aB2, sol.z));
alp3 = double(subs(alp3, sol.z));
end


function [rB, rD] = d_pos(phi1)
global AD BC rC, syms l_AB phi(t)
% set variables
rB = l_AB*[cos(phi(t)), sin(phi(t)), 0];
rD = AD*[cos(phi(t)), sin(phi(t)), 0];
% find rB
eq = (rB(1)-rC(1))^2 + (rB(2)-rC(2))^2 == BC^2;
sol = solve(eq, l_AB);
if double(subs(sol(1), phi, phi1))>0
    rB = subs(rB, l_AB, sol(1));
else
    rB = subs(rB, l_AB, sol(2)); end
end
function [vB, vD, omg3] = d_vel(phi1)
global rC
[rB, rD] = d_pos(phi1);
phi3 = atan((rB(2)-rC(2))/(rB(1)-rC(1)));
% find vB, vD, omg3
vB = diff(rB); vD = diff(rD);
omg3 = [0,0,diff(phi3)];
end
function [aB, aD, alp3] = d_acc(phi1)
[vB, vD, omg3] = d_vel(phi1);
% find aB, aD, alp3
aB = diff(vB); aD = diff(vD); alp3 = diff(omg3);
end


function [vB, vD, omg3, v12] = c_vel(phi1)
global rC omg1, syms x y z t
% set variables
[rB, rD] = pos(phi1);
omg23 = [0,0,x]; omg30 = [0,0,y]; v12 = [z,t,0];
% find omg23, omg30, a12
eq1 = omg1 + omg23 + omg30 == 0;
eq2 = cross(rB, omg23) + cross(rC, omg30) + v12 == 0;
eq3 = dot(cross(cross(rB, [1,0,0]), rB), v12) == 0;
sol = solve([eq1(3), eq2(1), eq2(2), eq3], [x,y,z,t]);
% return the results
omg23 = double(subs(omg23, x, sol.x));
omg30 = double(subs(omg30, y, sol.y));
v12 = double(subs(v12, [z,t], [sol.z,sol.t]));
omg3 = omg1 + omg23;
vB = cross(omg3, rB-rC);
vD = cross(omg1, rD);
end
function [aB, aD, alp3] = c_acc(phi1)
global rC omg1 alp1, syms x y z t
% set variables
[rB, rD] = pos(phi1); [vB, vD, omg2, v12] = c_vel(phi1);
alp23 = [0,0,x]; alp30 = [0,0,y]; a12 = [z,t,0];
% find alp23, alp30, a12
eq1 = alp1 + alp23 + alp30 == 0;
eq2 = cross(rB, alp23) + cross(rC, alp30) + a12 + 2*cross(omg1, v12) - omg2(3)^2*(rC-rB) - omg1(3)^2*rB == 0;
eq3 = dot(cross(cross(rB, [1,0,0]), rB), a12) == 0;
sol = solve([eq1(3), eq2(1), eq2(2), eq3], [x,y,z,t]);
% return the results
alp23 = double(subs(alp23, x, sol.x));
alp30 = double(subs(alp30, y, sol.y));
a12 = double(subs(a12, [z,t], [sol.z,sol.t]));
alp3 = alp1 + alp23;
aB = cross(alp3, rB-rC)-omg2(3)^2*(rB-rC);
aD = cross(alp1, rD)-omg1(3)^2*rD;
end
