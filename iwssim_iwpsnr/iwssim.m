function [iwssim iwmse iwpsnr]= iwssim(oimg, img, iw_flag, Nsc, K, L, weight, win, blk_size, parent, sigma_nsq)

MAX_PSNR = 1000;

if (nargin < 2 | nargin > 11)
    iwssim = -Inf;
    iwmse = -Inf;
    iwpsnr = -Inf;
    disp(['Input argument error']);
    return;
end

if (size(oimg) ~= size(img))
    iwssim = -Inf;
    iwmse = -Inf;
    iwpsnr = -Inf;
    disp(['Input argument error: images should have the same size']);
    return;
end

if (~exist('iw_flag'))
   iw_flag = 1;
end

if (~exist('Nsc'))
   Nsc = 5;
end

if (~exist('K'))
   K = [0.01 0.03]; % default from [Wang et al, IEEE-TIP 2004]
end

if (~exist('L'))
   L = 255;
end

if (~exist('weight'))
   weight = [0.0448 0.2856 0.3001 0.2363 0.1333]; % default from [Wang et al, IEEE-Asilomar 2003] 默认参数
end

if (exist('win')) 
    winSz = size(win); 
    if (winSz(1) ~= winSz(2)) % 保证滑动窗口是一个正方形
        iwssim = -Inf;
        iwmse = -Inf;
        iwpsnr = -Inf;
        disp(['Window size error']);
        return;
    end
    win_size = winSz(1);
else
    win_size = 11; % default from [Wang et al, IEEE-TIP 2004]
    gwin_std = 1.5; % default from [Wang et al, IEEE-TIP 2004]
    win = fspecial('gaussian', win_size, gwin_std); % default from [Wang et al, IEEE-TIP 2004]
    win = win/sum(win(:));
end

if (~exist('blk_size'))
   blk_size = 3; % spatial neighborhood size 空间领域大小
end

if (~exist('parent'))
   parent = 1; % include parent neighbor 包括父母的邻居
end

if (~exist('sigma_nsq'))
   sigma_nsq = 0.4; % default from [Sheikh & Bovik, IEEE-TIP 2006] 默认该论文的参数
end

blSzX = blk_size; 
blSzY = blk_size;
Nsc = min(5,Nsc); % 为5
weight = weight(1:Nsc);
weight = weight/sum(weight);%.^2是矩阵中的每个元素都求平方，^2是求矩阵（此矩阵为方阵）的平方。
if (min(size(oimg)) < win_size*(2^(Nsc-1))) % 判断原始图像要大于 88
    iwssim = -Inf;
    iwmse = -Inf;
    iwpsnr = -Inf;
    disp(['Image size too small for ' num2str(Nsc) ' scale IW-SSIM evaluation.']); % disp 输出
    return;
end
bound = ceil((win_size-1)/2); %ceil函数：朝正无穷大方向取整 向上取整
bound1 = bound - floor((blSzX-1)/2); % 向下取整

[pyro,pind]= buildLpyr(double(oimg),Nsc); %在矩阵(或向量)IM上构造一个拉普拉斯金字塔。 指定层数
[pyrd,pind]= buildLpyr(double(img),Nsc);%在矩阵(或向量)IM上构造一个拉普拉斯金字塔。
% 运用到了该文件夹下的另一个函数 
[cs_map l_map se_map] = scale_quality_maps(pyro,pyrd,pind,Nsc,K,L,win); % 得到了每层的ms-ssim系数 这里的层并不是简单的ms-ssim层 而是 金字塔中的层，保证了缩小的 时候的图片的 清晰度
if (iw_flag)
    iw_map = info_content_weight_map(pyro,pyrd,pind,Nsc,parent,blSzX,blSzY,sigma_nsq); % %计算尺度1到nssc -1的信息内容权重图
end

for s = 1:Nsc
    se = se_map{s}; % 计算同一层 两个子带的差异
    cs = cs_map{s};% 得到每层的c*s
    if (s==Nsc)
        cs = cs.*l_map; 
    end
    if (iw_flag)
        if (s<Nsc)
            iw = iw_map{s};
            iw = iw(bound1+1:end-bound1, bound1+1:end-bound1);
        else
            iw = ones(size(cs));
        end
        se = se(bound+1:end-bound, bound+1:end-bound);  % end是代表最后一个元素的位置      
        wmcs(s) = sum(sum(cs.*iw))/sum(sum(iw)); % 每层乘上对应的权重值 因为最后一层不包括 所以就只需要4层
        wmse(s) = sum(sum(se.*iw))/sum(sum(iw));
    else
        wmcs(s) = mean2(cs);
        wmse(s) = mean2(se);
    end
end

iwmse = prod(wmse(1:Nsc).^weight(1:Nsc)); % 信息权重的mse值
iwpsnr = min(10*log10(255^2/iwmse), MAX_PSNR);% psnr的信息权重
iwssim = prod(wmcs(1:Nsc).^weight(1:Nsc));% ssim的信息权重
