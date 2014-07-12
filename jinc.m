function[out] = jinc(x)
%
% jinc function
%
% J1(2*pi*x)/x
% divide by zero
%
% locate non-zero elements of x
mask=(x~=0);
%Intialise output with pi (value for x=0)
out =pi*ones(size(x));
%computer output values for all other x
out(mask)=besslj(1,2*pi*x(mask))./x(mask));
end
