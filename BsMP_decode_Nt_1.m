function u_hat = BsMP_decode_Nt_1(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr)

thre = 2;
ds = 4;
len_row = sum(F(1,:));
len_col = sum(F(:,1));
f = zeros(M, M, M, K*Nr);
Ivf = log(1/M*ones(K*Nr, M, V*Nt));
%Ivf = zeros(K*Nr, M, V*Nt);
Ifv = zeros(K*Nr, M, V*Nt);
Ifv_sum = zeros(1,M,V*Nt);
index_BsMP = zeros(K*Nr,ds,len_row);


for k = 1:K*Nr % resourses
    ind = find(F(k,:)==1); % non-zero elements, paths
    for m1 = 1:M
        for m2 = 1:M
            for m3 = 1:M
                f(m1,m2,m3,k) = -(1/N0)*abs(y(k)-(CB(k,m1,ind(1))*h(k,ind(1))+CB(k,m2,ind(2))*h(k,ind(2))+CB(k,m3,ind(3))*h(k,ind(3))  ))^2;
            end
        end
    end
end
 

%迭代开始
for mpa_iter = 1:thre
    
    %资源节点更新：Igv
    for k = 1:K*Nr
        ind = find(F(k,:)==1);
        for m1 = 1:M
            sIfv = zeros(M*M,1);
            for m2 = 1:M
                for m3 = 1:M
                    sIfv((m2-1)*M+m3,1) = f(m1,m2,m3,k) + Ivf(k,m2,ind(2)) + Ivf(k,m3,ind(3));
                end
            end
            Ifv(k,m1,ind(1)) = max(sIfv);
        end
        
        for m2 = 1:M
            sIfv = zeros(M*M, 1);
            for m1 = 1:M
                for m3 = 1:M
                    sIfv((m1-1)*M+m3,1) = f(m1,m2,m3,k) + Ivf(k,m1,ind(1)) + Ivf(k,m3,ind(3));
                end
            end
            Ifv(k,m2,ind(2)) = max(sIfv);
        end
        
        for m3 = 1:M
            sIfv = zeros(M*M, 1);
            for m1 = 1:M
                for m2 = 1:M
                    sIfv((m1-1)*M+m2,1) = f(m1,m2,m3,k) + Ivf(k,m1,ind(1)) + Ivf(k,m2,ind(2));
                end
            end
            Ifv(k,m3,ind(3)) = max(sIfv);
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


for mpa_iter = thre+1:iter_mpa_num
    
    for k = 1:K*Nr
        ind = find(F(k,:) == 1);
        for i = 1:len_row
            [~,I_index] = sort(Ivf(k,:,ind(i)),'descend');
            index_BsMP(k,:,i) = I_index(1:ds);
        end
    end
    
    
    %资源节点更新：Igv
    for k = 1:K*Nr
        ind = find(F(k,:)==1);
        for m1 = 1:M
            sIfv = zeros(ds*ds,1);
            for m2 = 1:ds
                for m3 = 1:ds
                    sIfv((m2-1)*ds+m3,1) = f(m1,index_BsMP(k,m2,2),index_BsMP(k,m3,3),k) + Ivf(k,index_BsMP(k,m2,2),ind(2)) + Ivf(k,index_BsMP(k,m3,3),ind(3));
                end
            end
            Ifv(k,m1,ind(1)) = max(sIfv);
        end
        
        
        
        for m2 = 1:M
            sIfv = zeros(ds*ds,1);
            for m1 = 1:ds
                for m3 = 1:ds
                    sIfv((m1-1)*ds+m3,1) = f(index_BsMP(k,m1,1),m2,index_BsMP(k,m3,3),k) + Ivf(k,index_BsMP(k,m1,1),ind(1)) + Ivf(k,index_BsMP(k,m3,3),ind(3));
                end
            end
            Ifv(k,m2,ind(2)) = max(sIfv);
        end
        
        for m3 = 1:M
            sIfv = zeros(ds*ds,1);
            for m1 = 1:ds
                for m2 = 1:ds
                    sIfv((m1-1)*ds+m2,1) = f(index_BsMP(k,m1,1),index_BsMP(k,m2,2),m3,k) + Ivf(k,index_BsMP(k,m1,1),ind(1)) + Ivf(k,index_BsMP(k,m2,2),ind(2));
                end
            end
            Ifv(k,m3,ind(3)) = max(sIfv);
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