%% Chebyshev Convergence
%% Introduction
% In this lab, we aim to investigate the convergence of the interpolations
% of three functions at the Chebyshev extreme points for each function as
% omega varies. The three functions, f(x), g(x), and h(x) are listed below:
% 
% f(x) = sin(omega*x+1)
%
% g(x) = sin(omega*x+1)/(4(x-0.1)^2+1)
%
% h(x) = sin(omega*x+1)*|(x-0.1)^3|
%
% Particularly, we first explore the length of preconvergent phase, and
% then evaluate the given 2 theorems based on corresponding functions. The
% 2 theorems are:
% 
% (1). If f is analytic in the closed Bernstein rho-ellipse, then Ck =
% O(rho^(-k-1)) and approximate errors O(rho^(-k)).
%
% (2). If f has a vth derivative of "bounded variation" in [1,-1], then Ck
% = O(k^(-v-1)) and approximate errors O(k^(-v)).
%
% Therefore, in order to examine these 2 theorems, it is the top priority
% to match each function with one of them. We will do so exactly after part
% one is well-handled.
% 
% n stands for the maximum degree of the Chebyshev polynomial of
% interpolation, and k is any arbitrary degree from 0 to n, inclusive:
%
n = 200;
k = 0:n;
%%
% Since we are interested in Chebyshev extreme points as the nodes, we have
% (pi * (0:n)) inside the cosine paratheses, where (0:n) is exactly k, the
% degree:
%
x = cos(pi * (0:n)' / n);
f = @(x,omega) sin(omega * x + 1);
g = @(x,omega) sin(omega * x + 1) ./ (4 * (x - 0.1).^ 2 + 1);
h = @(x,omega) sin(omega * x + 1) .* abs((x - 0.1).^ 3);
%%
% In part 1, we will use the three values in omega_values, as we'd like to
% see how variation of omega affects the length of preconvergent phase of
% each function. In part 2, we will use 20 as the omega value, as we intend
% to apply the theorems to our functions under a fixed omega, otherwise the
% result will be much less convincing:
omega_values = [100,60,20];
omega_value = 20;
%% PART I. Varying Omega and Length of Preconvergent Phase
% PART 1.1: Omega is 100
%
% Set up functions at omega = 100:
f1 = f(x,omega_values(1));
g1 = g(x,omega_values(1));
h1 = h(x,omega_values(1));
%%
% Set up the chebyshev coefficients of each function by dct1:
c_f1 = discrete_cosine_transform(f1)';
c_g1 = discrete_cosine_transform(g1)';
c_h1 = discrete_cosine_transform(h1)';
fig1 = figure;
semilogy(k,abs(c_f1), 'r', 'LineWidth',1, 'DisplayName', 'f(x)');
hold on;
semilogy(k,abs(c_g1), 'b', 'Linewidth',1, 'DisplayName', 'g(x)');
semilogy(k,abs(c_h1), 'g', 'LineWidth',1, 'DisplayName', 'h(x)');
xlabel('k, the degree');
ylabel('|Ck|, the abs. Cheb. Coefficients');
title(sprintf('Convergence Phase at \\omega = %d',omega_values(1)));
legend('f(x)','g(x)',"h(x)");
hold off;
%%
% Observation and formula: At omega = 100:
%
% The curves of all three functions maintain relatively stable at the
% beginning, and then decreases. f(x) and g(x) share pretty similar trend
% for which their curves, after the preconvergence phase, decreasew to some
% certain degrees, and then maintains relatively stable again during the
% roundoff plateau. However, in terms of h(x), although its curve also
% decreases after the preconvergence phase, it is overally smoother during
% the convergence phase, and we cannot observe the approximate moment when
% h(x) reaches the roundoff plateau.
%
% Formula: The preconvergent phase occurs for k=0 up to around omega.
%%
%
% PART 1.2: Omega is 60
f1 = f(x,omega_values(2));
g1 = g(x,omega_values(2));
h1 = h(x,omega_values(2));
c_f1 = discrete_cosine_transform(f1)';
c_g1 = discrete_cosine_transform(g1)';
c_h1 = discrete_cosine_transform(h1)';
fig2 = figure;
semilogy(k,abs(c_f1), 'r', 'LineWidth',1, 'DisplayName', 'f(x)');
hold on;
semilogy(k,abs(c_g1), 'b', 'Linewidth',1, 'DisplayName', 'g(x)');
semilogy(k,abs(c_h1), 'g', 'LineWidth',1, 'DisplayName', 'h(x)');
xlabel('k, the degree');
ylabel('|Ck|, the abs. Cheb. Coefficients');
title(sprintf('Convergence Phase at \\omega = %d',omega_values(2)));
legend('f(x)','g(x)',"h(x)");
hold off;
%%
% Observation and formula: At omega = 60:
%
% The key observation at omega=60 is similar to that at omega=100. The
% curves maintain relatively stable at first and then decrease. The curves
% of f(x) and g(x) drastically decreases during the convergence phase and
% finally exhibits relative stability again, whereas the curve of h(x)
% keeps decreasing as k increases up to n.
%
% Formula: The preconvergent phase occurs for k=0 up to around omega.
%%
%
% PART 1.3: Omega is 20
f1 = f(x,omega_values(3));
g1 = g(x,omega_values(3));
h1 = h(x,omega_values(3));
c_f1 = discrete_cosine_transform(f1)';
c_g1 = discrete_cosine_transform(g1)';
c_h1 = discrete_cosine_transform(h1)';
fig3 = figure;
semilogy(k,abs(c_f1), 'r', 'LineWidth',1, 'DisplayName', 'f(x)');
hold on;
semilogy(k,abs(c_g1), 'b', 'Linewidth',1, 'DisplayName', 'g(x)');
semilogy(k,abs(c_h1), 'g', 'LineWidth',1, 'DisplayName', 'h(x)');
xlabel('k, the degree');
ylabel('|Ck|, the abs. Cheb. Coefficients');
title(sprintf('Convergence Phase at \\omega = %d',omega_values(3)));
legend('f(x)','g(x)',"h(x)");
hold off;
%%
% Observation and formula: At omega = 20:
% 
% The key observation at omega=60 is similar to that at omega=100. The
% curves maintain relatively stable at first and then decrease. The curves
% of f(x) and g(x) drastically decreases during the convergence phase and
% finally exhibits relative stability again, whereas the curve of h(x)
% keeps decreasing as k increases up to n.
%
% Formula: The preconvergent phase occurs for k=0 up to around omega.
%% PART II. Convergence Phases and Bounding Curves
% Set up functions under fixed omega=20, and find the Chbyshev coefficients
% by dct1 what we did in part 1.
f2 = f(x,omega_value);
g2 = g(x,omega_value);
h2 = h(x,omega_value);
c_f2 = discrete_cosine_transform(f2)';
c_g2 = discrete_cosine_transform(g2)';
c_h2 = discrete_cosine_transform(h2)';
%%
% As we have noticed from part 1, when omega is 20, the degree where
% convergence phase starts is approximately 20, and therefore we set
% "kstart" as 20.
kstart = 20;
%%
% Part 2.1: f(x) and bounding curve at convergence phase
M_f = 10^11;
fig4 = figure;
semilogy(k,abs(c_f2), 'r', 'Linewidth',1, 'DisplayName', 'f(x)');
hold on;
semilogy((kstart:52), M_f ./ (3 .^(kstart:52)), 'b', 'Linewidth',1);
xlabel('k, the degree');
ylabel('|Ck|, the abs. Cheb. Coefficients');
title('f(x) Cheby. Coefficients against Degree and Bounding Curve');
legend('f(x) Cheby. Coefficients','Bounding Curve of f(x)')
hold off;
%%
% Analysis:
% 
% f(x) is sin(omega*x+1). Once we take advantage of its power series
% expansion and ratio test, we will see that the limit as n, the degree,
% goes to positive infinity, the ratio of a_(n+1)/a_(n) reaches 0. Ratio
% being 0 implies that f(x) is convergent for all x (infinite radius of 
% convergence, and thus analytic in the closed Bernstein rho-ellipse. 
% Therefore, we can take any arbitrary
% rho value as long as it appropriately fits the big-O approximation during
% convergence phase. Starting from 1.1, because rho must be greater than 1, 
% we eventually find that rho=3 is a suitable value of rho in this case.
%%
% Part 2.2: g(x) and bounding curve at convergence phase
M_g = 10^4;
z1 = ((0.2 + 1i) + sqrt(-4.96 + 0.4 * 1i)) / 2;
z2 = ((0.2 + 1i) - sqrt(-4.96 + 0.4 * 1i)) / 2;
z3 = ((0.2 - 1i) + sqrt(-4.96 - 0.4 * 1i)) / 2;
z4 = ((0.2 - 1i) - sqrt(-4.96 - 0.4 * 1i)) / 2;
rho =  abs(z1);
fig5 = figure;
semilogy(k,abs(c_g2), 'r', 'Linewidth',1, 'DisplayName', 'g(x)');
hold on;
semilogy((kstart:100), M_g ./ (rho .^(kstart:100)), 'b', 'Linewidth',1);
xlabel('k, the degree');
ylabel('|Ck|, the abs. Cheb. Coefficients');
title('g(x) Cheby. Coefficients against Degree and Bounding Curve');
legend('g(x) Cheby. Coefficients','Bounding Curve of g(x)')
hold off;
%%
% Analysis:
% 
% g(x) is sin(omega*x+1)/(4(x-0.1)^3+1). We notice that the numerator is
% convergent anywhere, and thus we pay attention to the denominator. If the
% denominator is 0, then g(x) will have a pole in the closed Bernstein
% rho-ellipse, although g(x) is still analytic since the denominator should 
% not be 0.Therefore, we let 4(x-0.1)^2+1=0, and by solving this 
% quadratic equation we get x=0.1+/-(i/2)=z0, where z0 denotes the pole.
% Next, we plug the x values into the equation: z0 = (z+(1/z))/2 and solve
% for z. We finally get 4 values of z as listed:
%
% z1 = ((0.2 + 1i) + sqrt(-4.96 + 0.4 * 1i)) / 2;
%
% z2 = ((0.2 + 1i) - sqrt(-4.96 + 0.4 * 1i)) / 2;
%
% z3 = ((0.2 - 1i) + sqrt(-4.96 - 0.4 * 1i)) / 2;
%
% z4 = ((0.2 - 1i) - sqrt(-4.96 - 0.4 * 1i)) / 2;
%
% Among these four, once we take their absolute values, |z1|=|z3|>1 and
% |z2|=|z4|<1. Because rho must be greater than 1, we choose |z1| as the 
% appropriate rho value.
%%
% Part 2.3: h(x) and bounding curve at convergence phase
M_h = 5;
fig6 = figure;
semilogy(k,abs(c_h2), 'r', 'LineWidth',1, 'DisplayName', 'h(x)');
hold on;
semilogy((kstart:180), M_h ./ ((kstart:180) .^2), 'b', 'LineWidth',1);
xlabel('k, the degree');
ylabel('|Ck|, the abs. Cheb. Coefficients');
title('h(x) Cheby. Coefficients against Degree and Bounding Curve');
legend('h(x) Cheby. Coefficients','Bounding Curve of h(x)')
hold off;
%%
% Analysis:
%
% h(x) is sin(omega*x+1)*|(x-0.1)^3|. As shown, the first part is exactly
% f(x), and so we pay attention to the second part as we did in part 2.2.
% Let p(x) denote the second absolute function.
% |(x-0.1)^3|=|x-0.1|^3 is an absolute function. Since absolute function
% has no general power series representation, we shall evaluate p(x)
% regarding to how many times it is differentiable, and thus determine the
% value for v. Since p(x) can be expressed as a piecewise function:
%
% p(x)=(x-0.1)^3, x>=0.1; p(x)=-(x-0.1)^3, x<0.1.
%
% First derivative:
%
% p'(x)=3(x-0.1)^2, x>0.1; p'(x)=-3(x-0.1)^2, x<0.1. Because both one-sided
% derivatives equal 0 at x=0.1, p'(0.1) exists.
%
% Second derivative:
%
% p''(x)=6(x-0.1), x>0.1; p''(x)=-6(x-0.1), x<0.1. Similarly, since both
% one-sided derivatives equal 0 at x=0.1, p'(0.1) exists.
%
% However, when we reach third derivative, the left-handed one is -6 when
% x<0.1, while the right-handed one is 6 when x>0.1. Since the two
% one-sided derivatives are different at x=0.1, p(x) is no longer
% differentiable. Therefore, p(x) is twice differentiable, and v is 2,
% which is the rho value for h(x).
%
% Noticeably from the graph, the bounding curve does not agree with the
% convergence phase of h(x) that well. It is a reasonable phenomenon since
% the bounding curve represents a big-O approximation, which is the upper
% bound. Besides, since the third derivative of p(x)=|(x-0.1)^3| has a jump
% at x=0.1, the discontinuity might impact the total variation of h(x) and
% thus make the bounding curve diverge from the convergence phase.