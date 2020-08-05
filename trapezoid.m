syms theta
num = 361;
ang = 6; b1 = 55; b2 = 5; b3 = 55; b4 = 245;
n = 500; omg = n*pi/30;
h = 2*cosd(ang);

s1 = trap(h, b1); s3 = -trap(h, b3) + h;
ds1 = diff(s1); ds3 = diff(s3);

n1 = fix(num*b1/360);
n2 = fix(num*b2/360);
n3 = fix(num*b3/360);
n4 = fix(num*b4/360);

s1n = zeros(n1,1); ds1n = zeros(n1,1);
s2n = zeros(n2,1); ds2n = zeros(n2,1);
s3n = zeros(n3,1); ds3n = zeros(n3,1);
s4n = zeros(n4,1); ds4n = zeros(n4,1);

len = n1 + n2 + n3 + n4;
T = linspace(0, 2*pi, len)';

for i=1:len
    if 1 <= i && i < n1 + 1
        s1n(i) = double(subs(s1, theta, T(i)));
        ds1n(i) = double(subs(ds1, theta, T(i)));
    elseif n1+1 <= i && i < n1+n2+1
        s2n(i-n1) = h;
    elseif n1+n2+1 <= i && i < n1+n2+n3+1
        s3n(i-n1-n2) = double(subs(s3, theta, T(i)-deg2rad(b1+b2)));
        ds3n(i-n1-n2) = double(subs(ds3, theta, T(i)-deg2rad(b1+b2)));
    end
end

sn = [s1n; s2n; s3n; s4n]; dsn = [ds1n; ds2n; ds3n; ds4n];
% plot(T, sn, 'linewidth', 1)
% plot(T, dsn, 'linewidth', 2)

s = sn; x = dsn; theta = T; Rb = 12;
RA = sqrt((Rb+s).^2 + x.^2);
sigA = theta + atan(x./(Rb + s));

xn = RA .* cos(sigA);
yn = -sign(omg) * RA .* sin(sigA);
axis equal
plot(xn,yn), pbaspect([1 1 1])


% function f = trap(h, beta)
%     syms t
%     beta = beta * pi / 180;
%     f = piecewise((0<=t & t<beta/8),...
%         h * (.38898448*t/beta - .0309544*sin(4*pi*t/beta)), ...
%         (beta/8<=t & t<3*beta/8), ...
%         h * (2.44406184*(t/beta)^2 - .22203097*t/beta + .00723407), ...
%         (3*beta/8<=t & t<5*beta/8),...
%         h * (1.6110154*t/beta - .0309544*sin(4*pi*t/beta-pi) - .3055077),... 
%         (5*beta/8<=t & t<7*beta/8),...
%         h * (-2.44406184*(t/beta)^2 + 4.6660917*t/beta - 1.2292648),...
%         (7*beta/8<=t & t<=beta),...
%         h * (.6110154 + .38898448*t/beta + .0309544*sin(4*pi*t/beta - 3*pi)));
% end

function f = trap(h, beta)
    syms theta
    beta = beta * pi / 180;
    f = piecewise((0<=theta & theta<beta/8),...
        0.09724613*h*(4*theta/beta-1/pi*sin(4*pi*theta/beta)), ...
        (beta/8<=theta & theta<3*beta/8), ...
        h*(2.444016188*(theta/beta)^2-0.22203094*theta/beta+0.00723406), ...
        (3*beta/8<=theta & theta<beta/2),...
        h*(1.6110154*theta/beta-0.0309544*sin(4*pi*theta/beta-pi)-0.3055077),... 
        (beta/2<=theta & theta<5*beta/8),...
        h*(1.6110155*theta/beta+0.0309544*sin(4*pi*theta/beta)-0.3055077),...
        (5*beta/8<=theta & theta<7*beta/8),...
        h*(4.6660917*theta/beta-2.44406188*(theta/beta)^2-1.2292650),...
        (7*beta/8<=theta & theta<=beta),...
        h*(0.6110155+0.3889845*theta/beta-0.0309544*sin(4*pi*theta/beta)));
end