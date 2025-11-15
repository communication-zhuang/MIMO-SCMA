function [P, Z] = Get_5GLDPC(n,R)

Zsize = [27, 54, 81];   % Square submatrices available size

if R == 1/2
    row = 12;           % Number of submatrices in a column
    rStr = '12';        % Rate in string form (12 = 1/2)
elseif R == 2/3
    row = 8;
    rStr = '23';
elseif R == 3/4
    row = 6;
    rStr = '34';
else
    row = 4;
    rStr = '56';
end



P = load(['protoH/',num2str(n),'_',rStr]);



if n == 648
    Z = Zsize(1);
elseif n == 1296
    Z = Zsize(2);
else
    Z = Zsize(3);
end





end