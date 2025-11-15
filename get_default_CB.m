function [CB,F,K,M,V] = get_default_CB(filetype)
switch filetype
    case 1
        load('Codebook\huawei.mat', 'CB');
    case 2
        load('Codebook\AMICB1.mat', 'CB');
    case 3
        load('Codebook\AMICB2.mat', 'CB');
end
 F = [0 1 1 0 1 0;
     1 0 1 0 0 1;
     0 1 0 1 0 1;
     1 0 0 1 1 0;];

 [K,M,V] = size(CB);
end


