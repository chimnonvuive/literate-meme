function out1 = yn(t)
%YN
%    OUT1 = YN(T)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    04-Aug-2020 12:06:50

if ((t < pi.*(1.0./2.88e2)) & (0.0 <= t))
    out1 = sin(t.*1.44e2).*(-1.54772e-1)+(t.*7.00172064e1)./pi;
elseif ((t < pi.*(1.0./9.6e1)) & (pi.*(1.0./2.88e2) <= t))
    out1 = (t.*(-3.99655746e1))./pi+t.^2.*1.0./pi.^2.*1.58375207232e4+3.617035e-2;
elseif ((t < pi.*(5.0./2.88e2)) & (pi.*(1.0./9.6e1) <= t))
    out1 = sin(t.*1.44e2).*1.54772e-1+(t.*2.89982772e2)./pi-1.5275385;
elseif ((t < pi.*(7.0./2.88e2)) & (pi.*(5.0./2.88e2) <= t))
    out1 = (t.*8.39896506e2)./pi-t.^2.*1.0./pi.^2.*1.58375207232e4-6.146324;
elseif ((t <= pi.*(1.0./3.6e1)) & (pi.*(7.0./2.88e2) <= t))
    out1 = sin(t.*1.44e2).*(-1.54772e-1)+(t.*7.00172064e1)./pi+3.055077;
else
    out1 = NaN;
end
