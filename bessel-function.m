%create custom Bessel function
lwr = 0;
uppr = pi;
theta = linspace(0,pi,N+1);

for li = 1:length(x)
% create the function to be intergrated
A = cos(x(ii)*sin(theta)-m*theta)/pi;

%^Integrate with Newton-coates method
J_m(ii) = nc_method(A,lwr,uppr,n);
end

%difference between matlab and custom
compare = J-J_m

v_eps = [eps, eps];
x_eps = [x(1),x(end)];

figure(1),
subplot(211),plot(x,J,x,J_m,':')
title('Comparision between MATLAB and custom bessel function')

legend('matlab', 'custom')
subplot(212),plot(x,compare,'-o',x_eps,v_eps,x_eps,-v_eps)
title('Difference between MATLAB and 
legend('diff' , '+eps', '-eps')
