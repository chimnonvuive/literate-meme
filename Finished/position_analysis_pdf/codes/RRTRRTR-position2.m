% find rF
BC_angle = atan((rC(2)-rD(2)) / (rC(1)-rD(1)));
rF = double(rC + pm*(DF-DC)*[cos(BC_angle), sin(BC_angle), 0]);

% find rG
D_angle = atan(rD(2) / rD(1));
rG = double(-pm*AG*[cos(D_angle), sin(D_angle), 0]);
end
