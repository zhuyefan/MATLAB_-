function imu = imenlarge2(im) % 将图片放大4倍后 再缩小2倍  变为原来的2倍

[M N] = size(im);
t1 = imresize(im, [4*M-3 4*N-3], 'bilinear'); %'bilinear'双线性插值 放大到指定大小
t2 = zeros(4*M-1, 4*N-1);
t2(2:end-1, 2:end-1) = t1; % 将t1 赋值过去t2
t2(1,:) = 2*t2(2,:) - t2(3,:);
t2(end,:) = 2*t2(end-1,:) - t2(end-2,:);
t2(:,1) = 2*t2(:,2) - t2(:,3);
t2(:,end) = 2*t2(:,end-1) - t2(:,end-2);
imu = t2(1:2:end, 1:2:end);  % 