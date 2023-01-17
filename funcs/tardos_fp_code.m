


% syms x
% int( 1 / (  sqrt(x) *  sqrt(1-x) )  , 0, x )
%
% -4*atan(   (   (1 - x)^(1/2) - 1    ) /x^(1/2)     )

function [Tardos_code,p]    = tardos_fp_code(epsilon, beta1, c0,c1)

L = ceil(4*pi^2*c0^2*log(1/beta1)); % c0 determines length

Tardos_code = zeros(c0,L);

dx = 10^8*eps;
x = 10^8*eps : dx : 1-10^7*eps; % limits depending on the random variable definition.

t =  mean (    [1/( exp(epsilon/2)+1 ) 1/c0]   );
f = PDF(x,t);
F = cumsum(f) * dx; % CDF
F = F / F(end); % recommended when max(F) is close to 1. Otherwise increase x - points to achieve F(end) - max value - close to 1.
% N = 30;
U = rand(1,c1); % c1 determines number of tardos code sequence
pdf2rand = @(u) x(find(u <= F, 1, 'first'));
p = arrayfun(pdf2rand, U);

for i = 1:c1
    Tardos_code(i,:) = rand(1,L)<p(i);
end

end

function f = PDF(x,t)
% Probability Distribution Function
f = 1/(  2*asin(1-2*t)  )  *   1 ./ (  sqrt(x) .*  sqrt(1-x) );
end