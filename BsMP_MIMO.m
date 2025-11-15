function u_hat = BsMP_MIMO(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr,index_MPA)

ds = 3;
len_row = sum(F(1,:));
len_col = sum(F(:,1));
f_k = zeros([M*ones(1,len_row) K*Nr]);
Ivf = log(1/M*ones(K*Nr, M, V*Nt));
Ifv = zeros(K*Nr, M, V*Nt);
Ivf_sum = zeros([M*ones(1,len_row) K*Nr]);
Ifv_sum = zeros(1,M,V*Nt);

for k = 1:K*Nr
    ind = find(F(k,:) == 1);
    h_temp = h(k,ind(1:len_row));
    for i = 1:M^(len_row)
        CB_temp = CB(index_MPA(k,:,i));
        sum_CB_h = sum(CB_temp.*h_temp,2);
        f_k(i + (k-1)*M^(len_row)) = -(1/N0)*abs(y(k) - sum_CB_h).^2;
    end   
end


for mpa_iter = 1:iter_mpa_num
    %资源节点更新：Ifv
    for k = 1:K*Nr            %MPA求和
        ind = find(F(k,:)==1);
        for i = 1:len_row
          [~,index_sort] = sort(Ivf(k,:,ind(i)),'descend');
        end
    end
    
    
    for k = 1:K*Nr            %MPA求和
        ind = find(F(k,:)==1);
        for i = 1:M^(len_row)
            Ivf_sum_temp = sum(Ivf(index_MPA(k,:,i)));
            Ivf_sum(i + (k-1)*M^(len_row)) = Ivf_sum_temp;
        end      
    end
    Ivf_f_k_sum = Ivf_sum + f_k;
  
    for k = 1:K*Nr            %求最大值
        ind = find(F(k,:)==1);
        for i = 1:len_row
            Ivf_f_k_sum_p = permute(reshape(Ivf_f_k_sum((k-1)*M^len_row+1:k*M^len_row),M*ones(1,len_row)),[i 1:i-1 i+1:len_row]);
            Ifv(k,:,ind(i)) =  max( (Ivf_f_k_sum_p - Ivf(k,:,ind(i))'),[],[2:len_row] )';
        end  
    end


    %变量节点更新：Ivg
    for k =1:V*Nt
        ind = find(F(:,k)==1);
        Ifv_sum(1,:,k) = sum( Ifv(ind(1:len_col),:,k) );
        % analogue of normalization in MPA, it can be removed (s1 and s2), but at high SNR and/or number of iterations NaN LLR values can be exist, so Max-Log-MPA is required
        for i = 1:len_col
            Ivf(ind(i),:,k) = Ifv_sum(1,:,k) - Ifv(ind(i),:,k);
            %Ivf(ind(i),:,k) = log( exp( Ivf( ind(i),:,k ) )./sum( exp( Ivf( ind(i),:,k ) ) ) );           
        end
    end
end

%最后的符号概率：Q
Q = zeros(M, V*Nt);
Q = squeeze(Ifv_sum(1,:,:));

% Q = log(exp(Q)./sum(exp(Q)));
u_hat = zeros(V*Nt,log2(M));
for k=1:V*Nt
    [~,L] = max(Q(:,k));
    u_hat(k,:) = Code(L,:);
end
end