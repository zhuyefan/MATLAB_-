function [psnr, mse] = psnr_mse(img1, img2, L)

MAX_PSNR = 1000;

if (nargin < 2 || nargin > 3) %nargin是用来判断输入变量个数的函数，特别是在利用了可变参数列表的函数中， 用nargin获取输入参数个数很方便。
   psnr = -Inf;
   mse = -Inf;
   disp('Error in psnr_mse.m: wrong input argument'); % 将这一行输出
   return;
end

if (size(img1) ~= size(img2))
   psnr = -Inf;
   mse = -Inf;
   disp('Error in psnr_mse.m: images need to have the same size');
   return;
end

if (~exist('L'))
   L = 255;
end

mse = mean2((double(img1) - double(img2)).^2); % .^2是矩阵中的每个元素都求平方  mean2 相当于mean( mean( A ) )相当于对整一个矩阵求像素平均值
psnr = 10*log10(L^2/mse);
psnr = min(MAX_PSNR, psnr);