global AB AC DC DF AG rA rC omg1 alp1
AB = .15; AC = .1; DC = .15; DF = .4; AG = .3;
rA = [0,0,0]; rC = [0,AC,0]; omg1 = [0,0,5/3*pi]; alp1 = [0,0,0];
% 
% hold on, axis equal
% [rB, rD, rF, rG] = pos(pi/6);
% plot([rA(1), rB(1)], [rA(2), rB(2)], 'r-o');
% plot([rD(1), rF(1)], [rD(2), rF(2)], 'g-o');
% plot([rA(1), rG(1)], [rA(2), rG(2)], 'b-o');
% 
% for phi1=0:2*pi/32:2*pi
%     hold on, axis equal
%     [rB, rD, rF, rG] = pos(phi1);
%     plot([rA(1), rB(1)], [rA(2), rB(2)], 'r-o');
%     plot([rD(1), rF(1)], [rD(2), rF(2)], 'g-o');
%     plot([rA(1), rG(1)], [rA(2), rG(2)], 'b-o');
% end
% 
% i=0; phis=zeros(201, 1);
% aGs = phis; vGs = phis;
% 
% for phi1=0:(2*pi/200):2*pi
%     i=i+1; phis(i)=phi1;
%     [~,~,~,vG,~,~] = vel(phi1); vGs(i,1) = vG(1);
%     [~,~,~,aG,~,~] = acc(phi1); aGs(i,1) = aG(1);
% end
% plot(phis, vGs(:,1))
% plot(phis, aGs(:,1))
% 
% phi1 = pi/6;
% syms t phi(t)
% old = {diff(phi, t, 2), diff(phi), phi}; new = {alp1(3), omg1(3), phi1};
% 
% [rB, rD, rF, rG] = pos(phi1)
% [d_rB, d_rD, d_rF, d_rG] = d_pos(phi1);
% d_rB = double(subs(d_rB, phi, phi1))
% d_rD = double(subs(d_rD, phi, phi1))
% d_rF = double(subs(d_rF, phi, phi1))
% d_rG = double(subs(d_rG, phi, phi1))
% 
% 
% [vB, vD, vF, vG, omg3, omg5] = vel(phi1)
% [c_vB, c_vD, c_vF, c_vG, c_omg3, c_omg5] = c_vel(phi1)
% [d_vB, d_vD, d_vF, d_vG, d_omg3, d_omg5] = d_vel(phi1);
% d_vB = double(subs(d_vB, old, new))
% d_vD = double(subs(d_vD, old, new))
% d_vF = double(subs(d_vF, old, new))
% d_vG = double(subs(d_vG, old, new))
% d_omg3 = double(subs(d_omg3, old, new))
% d_omg5 = double(subs(d_omg5, old, new))
% 
% [aB, aD, aF, aG, alp3, alp5] = acc(phi1)
% [c_aB, c_aD, c_aF, c_aG, c_alp3, c_alp5] = c_acc(phi1)
% [d_aB, d_aD, d_aF, d_aG, d_alp3, d_alp5] = d_acc(phi1);
% d_aB = double(subs(d_aB, old, new))
% d_aD = double(subs(d_aD, old, new))
% d_aF = double(subs(d_aF, old, new))
% d_aG = double(subs(d_aG, old, new))
% d_alp3 = double(subs(d_alp3, old, new))
% d_alp5 = double(subs(d_alp5, old, new))


function [rB, rD, rF, rG] = pos(phi1)
global AB DC DF AG rC
% set variables
rB = AB*[cos(phi1), sin(phi1), 0];
pm = sign(rB(1));
% find rD
CB_angle = atan((rB(2)-rC(2)) / (rB(1)-rC(1)))+pi;
rD = double(rC + pm*DC*[cos(CB_angle), sin(CB_angle), 0]);
% find rF
BC_angle = atan((rC(2)-rD(2)) / (rC(1)-rD(1)));
rF = double(rC + pm*(DF-DC)*[cos(BC_angle), sin(BC_angle), 0]);
% find rG
D_angle = atan(rD(2) / rD(1));
rG = double(-pm*AG*[cos(D_angle), sin(D_angle), 0]);
end
function [vB, vD, vF, vG, omg3, omg5, vB3B2, vD5D4] = vel(phi1)
global omg1 rC, syms x y z
% set variables
[rB, rD, rF, rG] = pos(phi1);
omg3 = [0,0,z]; vB3B2 = [x,y,0];
omg5 = [0,0,z]; vD5D4 = [x,y,0];
vB = cross(omg1, rB);
vD = cross(omg3, rD-rC);
vF = cross(omg3, rF-rC);
vG = cross(omg5, rG);
% find  vB3B2, omg3, vD, vF
eq1 = cross(omg3, rB-rC) == vB + vB3B2;
eq2 = dot(cross(cross(rC-rB, [1,0,0]), rC-rB), vB3B2) == 0;
sol = solve([eq1(1),eq1(2),eq2], [x,y,z]);
% return the values
vB3B2 = double(subs(vB3B2, [x, y], [sol.x, sol.y]));
omg3 = double(subs(omg3, sol.z));
vD = double(subs(vD, sol.z));
vF = double(subs(vF, sol.z));
% find vB5B4, omg5, vG
eq1 = cross(omg5, rD) == vD + vD5D4;
eq2 = dot(cross(cross(rG, [1,0,0]), rG), vD5D4) == 0;
sol = solve([eq1(1),eq1(2),eq2], [x,y,z]);
% return the values
vD5D4 = double(subs(vD5D4, [x, y], [sol.x, sol.y]));
omg5 = double(subs(omg5, sol.z));
vG = double(subs(vG, sol.z));
end
function [aB, aD, aF, aG, alp3, alp5] = acc(phi1)
global rC omg1 alp1, syms x y z
% set variables
[rB, rD, rF, rG] = pos(phi1);
[~,~,~,~, omg3, omg5, vB3B2, vD5D4] = vel(phi1);
alp3 = [0,0,x]; aB3B2r = [y,z,0];
alp5 = [0,0,x]; aD5D4r = [y,z,0];
aB = cross(alp1, rB) - norm(omg1)^2*rB;
aD = cross(alp3, rD-rC) - norm(omg3)^2*(rD-rC);
aF = cross(alp3, rF-rC) - norm(omg3)^2*(rF-rC);
aG = cross(alp5, rG) - norm(omg5)^2*rG;
aB3B2c = 2*cross(omg3, vB3B2);
aD5D4c = 2*cross(omg5, vD5D4);
% find alp3, aB3B2r
eq1 = cross(alp3, rB-rC) - norm(omg3)^2*(rB-rC) == aB + aB3B2c + aB3B2r;
eq2 = dot(cross(cross(rC-rB, [1,0,0]), rC-rB), aB3B2r) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);
% return the values
alp3 = double(subs(alp3, sol.x));
aB3B2r = double(subs(aB3B2r, [y,z], [sol.y,sol.z]));
aD = double(subs(aD, sol(1)));
aF = double(subs(aF, sol(1)));
% find alp5, aD5D4r
eq1 = cross(alp5, rD) - norm(omg5)^2*(rD) == aD + aD5D4c + aD5D4r;
eq2 = dot(cross(cross(rG, [1,0,0]), rG), aD5D4r) == 0;
sol = solve([eq1(1), eq1(2), eq2], [x,y,z]);
% return the values
alp5 = double(subs(alp5, sol.x));
aD5D4r = double(subs(aD5D4r, [y,z], [sol.y,sol.z]));
aG = double(subs(aD, sol(1)));
end

function [rB, rD, rF, rG] = d_pos(phi1)
global AB DC DF AG rC, syms phi(t) x y
rB = AB*[cos(phi(t)), sin(phi(t)), 0];
pm = sign(subs(rB(1), phi, phi1));
% find rD
CB_angle = atan((rC(2)-rB(2))/(rC(1)-rB(1))) + pi;
rD = rC + pm*DC*[cos(CB_angle), sin(CB_angle), 0];
% find rF
DC_angle = atan((rD(2)-rC(2)) / (rD(1)-rC(1)));
rF = rC + pm*(DF-DC)*[cos(DC_angle), sin(DC_angle), 0];
% find rG
D_angle = atan(rD(2)/rD(1));
rG = -pm*AG*[cos(D_angle), sin(D_angle), 0];
end
function [vB, vD, vF, vG, omg3, omg5] = d_vel(phi1)
global rC
% set variables
[rB, rD, rF, rG] = d_pos(phi1);
phi3 = atan((rD(2)-rC(2))/(rD(1)-rC(1)));
phi5 = atan(rG(2)/rG(1));
% find vB, vD, vF, vG, omg3, omg5
vB = diff(rB); vD = diff(rD); vF = diff(rF); vG = diff(rG);
omg3 = [0,0,diff(phi3)]; omg5 = [0,0,diff(phi5)];
end
function [aB, aD, aF, aG, alp3, alp5] = d_acc(phi1)
[vB, vD, vF, vG, omg3, omg5] = d_vel(phi1);
aB = diff(vB); aD = diff(vD); aF = diff(vF); aG = diff(vG);
alp3 = diff(omg3); alp5 = diff(omg5);
end


function [vB, vD, vF, vG, omg3, omg5, v23, v45] = c_vel(phi1)
global rC omg1, syms x y z t
% set variables
[rB, rD, rF, rG] = pos(phi1); vB = cross(omg1, rB);
omg12 = [0,0,x]; v23 = [y,z,0]; omg30 = [0,0,t];
omg34 = [0,0,x]; v45 = [y,z,0]; omg50 = [0,0,t];
% find v23, omg3, vD, vF of contour #1
eq1 = omg1 + omg12 + omg30 == 0;
eq2 = dot(cross(cross(rC-rB, [1,0,0]), rC-rB), v23) == 0;
eq3 = cross(rB, omg12) + cross(rC, omg30) + v23 == 0;
sol = solve([eq1(3), eq2, eq3(1), eq3(2)], [x, y, z, t]);
% return the values
omg12 = double(subs(omg12, sol.x));
v23 = double(subs(v23, [y,z], [sol.y,sol.z]));
omg30 = double(subs(omg30, sol.t));
omg3 = -omg30;
vD = cross(omg3, rD-rC);
vF = cross(omg3, rF-rC);
% find v45, omg5, vG of contour #2
eq1 = -omg30 + omg34 + omg50 == 0;
eq2 = dot(cross(cross(rG, [1,0,0]), rG), v45) == 0;
eq3 = cross(rC, -omg30) + cross(rD, omg34) + v45 == 0;
sol = solve([eq1(3),eq2,eq3(1),eq3(2)], [x,y,z,t]);
% return the values
omg34 = double(subs(omg34, sol.x));
v45 = double(subs(v45, [y, z], [sol.y, sol.z]));
omg50 = double(subs(omg50, sol.t));
omg5 = -omg50;
vG = cross(omg5, rG);
end
function [aB, aD, aF, aG, alp3, alp5] = c_acc(phi1)
global rC omg1 alp1, syms x y z t
% set variables
[rB, rD, rF, rG] = pos(phi1);
[~,~,~,~, omg3, omg5, v23, v45] = c_vel(phi1);
alp12 = [0,0,x]; a23 = [y,z,0]; alp30 = [0,0,t];
alp34 = [0,0,x]; a45 = [y,z,0]; alp50 = [0,0,t];
% find alp3, aD, aF of contour 1
eq1 = alp1 + alp12 + alp30 == 0;
eq2 = dot(cross(cross(rC-rB, [1,0,0]), rC-rB), a23) == 0;
eq3 = cross(rB, alp12) + cross(rC, alp30) + 2*cross(omg3, v23) - omg1(3)^2*rB - omg3(3)^2*(rB-rC) + a23 == 0;
sol = solve([eq1(3),eq2,eq3(1),eq3(2)], [x,y,z,t]);
% return the values
alp12 = double(subs(alp12, sol.x));
a23 = double(subs(a23, [y,z], [sol.y,sol.z]));
alp30 = double(subs(alp30, sol.t));
alp3 = -alp30;
aB = cross(alp1, rB) - omg1(3)^2*rB;
aD = cross(alp3, rD-rC) - omg3(3)^2*(rD-rC);
aF = cross(alp3, rF-rC) - omg3(3)^2*(rF-rC);
% find alp5, aG of contour 2
eq1 = -alp30 + alp34 + alp50 == 0;
eq2 = dot(cross(cross(rG, [1,0,0]), rG), a45) == 0;
eq3 = -cross(rC, alp30) + cross(rD, alp34) + 2*cross(omg5, v45) - omg3(3)^2*(rD-rC) - omg5(3)^2*rD + a45 == 0;
sol = solve([eq1(3),eq2,eq3(1),eq3(2)], [x,y,z,t]);
% return the values
alp34 = double(subs(alp34, sol.x));
a45 = double(subs(a45, [y,z], [sol.y,sol.z]));
alp50 = double(subs(alp50, sol.t));
alp5 = -alp50;
aG = cross(alp5, rG) - omg5(3)^2*rG;
end
