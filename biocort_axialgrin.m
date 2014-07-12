% AXIAL GRADIENT INDEX LENS ABERRATIONS SCRIPT
function [sph, coma, astig, curv, dist, Tch, Lch] = biocort(ie, d, z)

%sum of gradient index
ie = 10;
d=1;
h=1;
n0 = 
%EQ 9
%integrate
% ray propagation in an inhomogeneous media
% see Wikipedia gradient index optics
z = @(n) exp(-n^2);
%dzdivnz = diff(z);
dzdivnz = 1;
%int_d = sum(dzdivnz,d,0);
int_d=1;
n_1 = 1/d*int_d;

%EQ 10 
%int_d3 = integrate(dzdivnz.^3,d,0);
int_d3 = 1;
n_3 = 1/d*int_d3;

%EQ 14
%Transfer equation for marginal ray
h_dash = h - d*n_1*n0*u;
%EQ 15
n0_dash_u_dash=n0*u;

%calculate the primrary aberration coefficients
% p=1,2,3,4 are spherical, coma, astigmatism and distortion
lambda_p = sum(Sp + Sp_star) + sum(Tp);
% P = petzval curvature Ps(surfaces), Pt(GRIN media)
P = sum(Ps) + sum(Pt);

%chromatic aberrations for axial and lateral colour
Lambda_lambdap = sum(S_lambdap) + sum(T_lambdap);

% paraxial invariants homogeneous case
for i=1:ie
S1 = (n0*i)^2*h*((u(i)/n0(i))-(u(i-1)/n0(i-1)));
S2 = n0*i*n0*j*h*(u(i)/n0(i)-u(i-1)/n0(i-1));
S3 = (n0*j)^2*h*(u(i)/n0(i))-(u(i-1)/n0(i-1));
P4 = -rho*H^2;
S4 = (n0*j)^2*m*(u(i)-n0(i))-u(i-1)/n0(i-1);
S_l1 = h*n0*i*(dell_lam(i)*n0(i)/n0(i)-dell_lam(i-1)*(n0(i-1)/n0(i-1)));
S_lam2 = h*n0*j*((delta_lam(i)*n0(i)/n0(i)-delta_lam(i-1)*n0(i-1)/n0(i-1)));
end
delta_lam_n0 = n0_F-n_oC;

%including paraxial invariants
n0_i = n0*h*rho;
n0_j = n0*m*rho-n0*w;
H = m*n0*u-h*n0*w;

%EQ22
%first refractive index at surface vertex
%differentiate
N_z = dn/dz;

%inhomogeneou contribtuions to seidel coefficients are
for i = 1:ie
S1_star = h^4*rho^2*N_z(i)-N_z(i-1);
S2_star = h^3*m*rho^2*N_z(i)-N_z(i-1);
S3_star = h^2*m^2*rho^2*N_z(i)-N_z(i-1);
S4star = h*m^3*rho^2*N_z(i)-N_z(i-1);

%transfer contributions
T1 = n0^3*u^3*(n0*u*d*n_3+(h(i)/n0(i)^2)-(h(i-1)/n0(i-1)^2));
end
T4divT3 = (n0*w)/(n0*u);
T3divT2 = (n0*w)/(n0*u);
T2divT1 = (n0*w)/(n0*u);
Pt = 0;

%transfer contributions to the chromatic coefficients are

T_lam1 = n0*u*(n0*u*dell_lam*n_1-(dell_lam(i)*n0(i)/n0(i)*h(i)))-(dell_lam(i-1)*n0(i-1))/n0(i-1)*h(i-1);
T_lam2 = n0*w*(n0*u*d*dell_lam*n_1-(dell_lam(i)*n0(i)*h(i)/n0(i)-dell_lam(i-1)*n0(i-1)*h(i-1)/n0(i-1)));
dell_lam_n_1 = 1/nF - 1/nC;


%Save 'Answer Axial contributions to gradient index lens'
