[a,str] = xlsread('1.xlsx','A2:A15')
[a,str1] = xlsread('1.xlsx','B62:B75')
result=[1;1;1;1;1;1;1;1;1;1;1;1;1;1]
for i = 1:14
    x = imread(str{i});
    x = rgb2gray(x);
    y = imread(str1{i});
    y = rgb2gray(y);
    
    result(i) =iwssim(x,y)
end

    

    
