function [rB, rC] = pos(phi1)
global AB BC
syms x

rB = AB*[cos(phi1), sin(phi1), 0];
rC = [AB*cos(phi1)+BC*asin(AB*sin(phi1)/BC), 0, 0];
end