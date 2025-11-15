function [c_hat,c_llr] = bp_LDPC_decode(llr,H,V,n,m,iter_LDPC)

itr = 0;
success = 0;
llr_0 = permute(llr,[3 1 2]);
llr_P = repmat(llr_0,m,1,1);
L_v2c = llr_P;
L_c2v = zeros(m,n,V);
sum_Lc2v = zeros(1,n,V);

while (itr<iter_LDPC) && (success == 0)
    itr = itr + 1;

    for i = 1:m
        idx_i = find(H(i,:));
        prod_Lv2c = prod( tanh(L_v2c(i,idx_i,:)/2),2 )./ tanh(L_v2c(i,idx_i,:)/2);
        prod_Lv2c = min(max(prod_Lv2c,-1 + 1e-16), 1-1e-16);
        L_c2v(i,idx_i,:) = 2*atanh(prod_Lv2c);

    end

    for j = 1:n
        idx_j = find(H(:,j));
        sum_Lc2v(1,j,:) = sum(L_c2v(idx_j,j,:),1);
        L_v2c(idx_j,j,:) = sum_Lc2v(1,j,:) - L_c2v(idx_j,j,:) + llr_P(idx_j,j,:);
    end

    c_llr = squeeze(sum_Lc2v) + llr;
    c_hat = (1 - sign(c_llr))/2;

    check = sum(mod(H*c_hat,2),'all');

    if check==0
        success = 1;
    end
   
end

end




