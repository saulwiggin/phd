%fraun_circ - Fraunhofer irradiance plot

L=0.2;
M=250;
dx=L/M;
x=-L/2:dx:L/2-dx; y=x;
[X,Y]=meshgrid(x,y);

w=1e-3;
lambda=0.633e-6;
z=50;
k=2*pi/lambda;
lz=lambda*z;

%irradiance
I2=(w^2/lz)^2.*(jinc(w/lz*sqrt(X.^2+Y^2))).^2;

figure(1)
imagesc(x,y,nthroot(I2,3));
xlabel('x (m)'); ylabel('y (m)');
colormap('gray');
axis square;
axis xy;

figure(2)
plot(x,I2(M/2+1,:));
xlabel('x(m)'); ylabel('Irradiance');
