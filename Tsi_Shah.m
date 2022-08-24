function z = Tsi_Shah(E,tilt, slant, niter)

% Ping-Sing Tsai and Mubarak Shah surface depth reconstruction method:
% 
% P.S. Tsai and M. Shah, "Shape from shading using linear approximation",
% Image and Vision Computing Journal, Vol. 12, No.8, pp. 487-498, 1994.
% 
% Inputs:
%  E - shading image
%  niter - the number of iterations
%  tilt, slant - illumination direction angles
% Output:
%  z - estimated depth

ps = cos(tilt)*tan(slant);
qs = sin(tilt)*tan(slant);

opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);

R = @(p,q) (cos(slant) + p .* cos(tilt)*sin(slant)+ q .* ...
        sin(tilt)*sin(slant))./sqrt(1 + p.^2 + q.^2);

z = zeros(size(E));
for n = 1:niter
    
    del = opD(z);
    p = del(:,:,2);
    q = del(:,:,1);
    Rout = R(p,q);
    
    Rout = max(0,Rout);
    f = E - Rout;
    
    df = (p + q).*(ps*p + qs*q + 1)./(sqrt((1 + p.^2 + q.^2).^3)* ...
        sqrt(1 + ps^2 + qs^2)) - (ps + qs)./(sqrt(1 + p.^2 + q.^2)* ...
        sqrt(1 + ps^2 + qs^2));
    
    z = z - f./(df + eps);
    
end
