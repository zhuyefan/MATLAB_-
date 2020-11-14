%contrast-structure similarity map and squared error map for each scale, and luminance similarity map for the coarsest scale
function [cs_map l_map se_map]= scale_quality_maps(pyro,pyrd,pind,Nsc,K,L,win)
    if (nargin < 3 | nargin > 7)
        cs_map = -Inf;
        l_map = -Inf;
    	se_map = -Inf;
        disp(['Input argument error']);
        return;
    end
    if (~exist('Nsc'))
        Nsc = size(pind,1);
    end
    if (~exist('K'))
        K = [0.01 0.03]; % default from [Wang et al, IEEE-TIP 2004]
    end
    if (~exist('L'))
        L = 255;
    end
    if (~exist('win'))
        win = fspecial('gaussian', 11, 1.5); % default from [Wang et al, IEEE-TIP 2004]
        win = win/sum(win(:));
    end

    pyro = real(pyro); % 复数的实部数值 得到拉普拉斯金字塔的高度精确值即可
    pyrd = real(pyrd);
    C1 = (K(1)*L)^2;
    C2 = (K(2)*L)^2;

    for i=1:Nsc
        y = pyrBand(pyro, pind, i); %原图  按顺序访问金字塔的一个子带(高斯，拉普拉斯算子，QMF/小波，或可操纵的)
        yn = pyrBand(pyrd, pind, i);% 失真图像 按顺序访问金字塔的一个子带(高斯，拉普拉斯算子，QMF/小波，或可操纵的)
        se_map{i} = (y-yn).^2; % 计算同一层 两个子带的差异
        
        mu1 = filter2(win, y, 'valid');
        mu2 = filter2(win, yn, 'valid');% 得到一个高斯滤波处理后的图像
        sigma12 = filter2(win, y.*yn, 'valid') - mu1.*mu2; % 协方差
        sigma1_sq  = filter2(win, y.^2, 'valid') - mu1.^2; % 方差
        sigma2_sq  = filter2(win, yn.^2, 'valid') - mu2.^2;% 方差
        sigma1_sq = max(0, sigma1_sq);
        sigma2_sq = max(0, sigma2_sq);
        cs_map{i} = (2*sigma12 + C2)./(sigma1_sq + sigma2_sq + C2); % 得到每层的c*s
        
        if (i == Nsc)
            l_map = (2*mu1.*mu2 + C1)./(mu1.^2 + mu2.^2 + C1); % 最后一层乘以l
        end
    end