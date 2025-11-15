function H = GetLDPC_H(n, R)

Zsize = [27, 54, 81];   % Square submatrices available size

col = 24;               % Number of submatrices in a row

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

% This block will be needed only if the support for different
% codeword size will be implemented
if n == 648
    Z = Zsize(1);
elseif n == 1296
    Z = Zsize(2);
else
    Z = Zsize(3);
end

H = zeros(row*Z,col*Z);

protoH = load(['protoH/',num2str(n),'_',rStr]);

for i = 1:row
    for j = 1:col
        if protoH(i,j) >= 0
            A = eye(Z);
            k = protoH(i,j);
            H((i-1)*Z+1:(i-1)*Z+Z,(j-1)*Z+1:(j-1)*Z+Z) = circshift(A,[0 k]);
        end
    end
end

end