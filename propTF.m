function[u2]=propTF(u1,L,lambda,z);
%propagation - transfer function approach
%assumes same x and y side lenegths and 
% uniform sampling
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field

[M,N] = size(u1);
dx=L/M;
k=2*pi/lambda;

fx=-1/(2*dx):1/L:1/(2*dx)-1/L;
[FX,FY]=meshgrid(f,fy);

Hexp(-j*pi*lambda*z*(FX.^2+FY.^2));
H=fftshift(H);
U1 = fft2(fftshift(u1));
U2 = H.*U1;
u2=ifftshift(ifft2(U2));
end
