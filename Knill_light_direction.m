

function [tilt, slant] = Knill_light_direction(E)

% Compute light direction tilt and slant angle using David C. Knill's
% method:
%
% Knill, D. C. (1990) Estimating illuminant direction and degree of surface relief. 
% Journal of the Optical Society of America A, 7 (4), 759-775
%
% Input:
%  - E : is a shading image

global E0;
E0 = E;

x0 = ones(1,2);
cost = @(x) lz_nz_C(x) + lz_nz_R(x);
[x,fval,exitflag,output] = fminsearch(cost,x0);

slant = acos(x(1));

D = [-0.0577 0.215 -0.804 0 0.804 - 0.215 0.0577];
dIdx = imfilter(E,D);
dIdy = imfilter(E,D');
tilt = 0.5*atan( (2*mean(dIdx(:).*dIdy(:))) / (mean(dIdx(:).^2) - mean(dIdy(:).^2)) );

function f = lz_nz_C(x)

global E0;

c0 = (mean(E0(:).^2) - mean(E0(:))^2)/mean(E0(:))^2;
lz = x(1);
sigma = x(2);

Enz = pi/(sqrt(2)*sigma)*exp(1/(2*sigma^2))*erfc(1/(sqrt(2)*sigma));
Enz2 = exp(1/(2*sigma^2))*expint(1/(2*sigma^2))/(2*sigma^2);
c = (1 - lz^2 + (3*lz^2 - 1)*Enz2)/(2*lz^2*Enz^2) - 1;

f = (c - c0)^2;

function f = lz_nz_R(x)

global E0

D = [-0.0577 0.215 -0.804 0 0.804 - 0.215 0.0577];
dIdx = imfilter(E0,D);
dIdy = imfilter(E0,D');

r0 = mean(dIdx(:).^2)/mean(dIdy(:).^2);

lz = x(1);
sigma = x(2);

Enz2 = exp(1/(2*sigma^2))*expint(1/(2*sigma^2))/(2*sigma^2);
Enz4 = (1 - Enz2)/(2*sigma^2);
Enz6 = (1 - Enz4)/(4*sigma^2);

r = (5*Enz2 + 2*Enz4 + 5*Enz6 - lz^2*(5*Enz2 - 6*Enz4 + 13*Enz6))/...
    (3*Enz2 - 2*Enz4 + 3*Enz6 - lz^2*(3*Enz2 - 10*Enz4 + 11*Enz6));

f = (r - r0)^2;
