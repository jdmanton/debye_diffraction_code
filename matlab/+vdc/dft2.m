function Xhat = dft2(X, k, a, N)
% 2D discret fourier transform
%
% Xhat = dft(X, k, a, N)
% Xhat = dft(X, k, a)
% X : 2 image
% k : frequency offsets  [row, columns]
% a : upsampling factor
% N : number of samples [row, columns]
%
% Jerome Boulanger & James Manton (2018)
%
if nargin < 4
    N = size(X);
end
if size(N) == 1
    N = [N N];
end
f1 = linspace(-1/(2*a(1)),1/(2*a(1)),N(1)) - k(1) / size(X,1);
f2 = linspace(-1/(2*a(2)),1/(2*a(2)),N(2)) - k(2) / size(X,2);
x1 = 0:size(X,1)-1; 
x2 = 0:size(X,2)-1;
F1 = exp(-1i*2*pi*f1'*x1);
F2 = exp(-1i*2*pi*x2'*f2);
Xhat = F1 * X * F2;
