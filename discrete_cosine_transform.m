%=======================================
% Discrete cosine transform
%---------------------------------------
% Code adapted from
% P. Moin, Fundamentals of Engineering Numerical Analysis (2010)
%=======================================
function y=discrete_cosine_transform(f_vals)

N      =length(f_vals);
y      =[f_vals;flipud(f_vals(2:N-1,:))];
y_ft   =fft(y);
y      =real(y_ft(1:N,:)/(N-1));
y(1,:) = y(1,:)/2;
y(N,:) = y(N,:)/2;
