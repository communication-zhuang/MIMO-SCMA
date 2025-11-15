function [u_llr,u_hat] = Max_Log_Nt_1_coded_vector(y,K,V,M,F,CB,N0,h,iter_mpa_num,Nt,Nr,slot,bit_num)
len_row = sum(F(1,:));
len_col = sum(F(:,1));
f = zeros(M, M, M, K*Nr,slot);
Ivf = log(1/M*ones(K*Nr, M, V*Nt,slot));
Ifv = zeros(K*Nr, M, V*Nt,slot);
Ifv_sum = zeros(1,M,V*Nt,slot);

for k = 1:K*Nr % resourses
    ind = find(F(k,:)==1); % non-zero elements, paths
    for m1 = 1:M
        for m2 = 1:M
            for m3 = 1:M
                f(m1,m2,m3,k,:) = -(1/N0)*abs(y(k,:)-(CB(k,m1,ind(1))*h(k,ind(1))+CB(k,m2,ind(2))*h(k,ind(2))+CB(k,m3,ind(3))*h(k,ind(3))  )).^2;
            end
        end
    end
end
 

%迭代开始
for mpa_iter = 1:iter_mpa_num
    Ivf0 = permute(Ivf,[1 2 3 5 4]);
    %资源节点更新：Igv
    for k0 = 1:K*Nr
        k = k0;
        ind = find(F(k,:)==1);
        for m1 = 1:M
            sIfv = zeros(M*M,1,slot);
            for m2 = 1:M
                for m3 = 1:M
                    sIfv((m2-1)*M+m3,1,:) = f(m1,m2,m3,k,:) + Ivf0(k,m2,ind(2),1,:) + Ivf0(k,m3,ind(3),1,:);
                    %sIfv((m2-1)*M+m3,1,:) = reshape(f(m1,m2,m3,k,:),[],slot) + reshape(Ivf(k,m2,ind(2),:),[],slot) + reshape(Ivf(k,m3,ind(3),:),[],slot);
                end
            end
            Ifv(k,m1,ind(1),:) = max(sIfv,[],1);
        end
        %Ifv(k,:,ind(1),:) = Ifv(k,:,ind(1),:) - Ifv(k,1,ind(1),:);
        
        for m2 = 1:M
            sIfv = zeros(M*M, 1,slot);
            for m1 = 1:M
                for m3 = 1:M
                    sIfv((m1-1)*M+m3,1,:) = f(m1,m2,m3,k,:) + Ivf0(k,m1,ind(1),1,:) + Ivf0(k,m3,ind(3),1,:);
                    %sIfv((m1-1)*M+m3,1,:) = reshape(f(m1,m2,m3,k,:),[],slot) + reshape(Ivf(k,m1,ind(1),:),[],slot) + reshape(Ivf(k,m3,ind(3),:),[],slot);
                end
            end
            Ifv(k,m2,ind(2),:) = max(sIfv,[],1);
        end
        %Ifv(k,:,ind(2),:) =  Ifv(k,:,ind(2),:) -  Ifv(k,1,ind(2),:);
        
        for m3 = 1:M
            sIfv = zeros(M*M, 1,slot);
            for m1 = 1:M
                for m2 = 1:M
                    sIfv((m1-1)*M+m2,1,:) = f(m1,m2,m3,k,:) + Ivf0(k,m1,ind(1),1,:) + Ivf0(k,m2,ind(2),1,:);
                    %sIfv((m1-1)*M+m2,1,:) = reshape(f(m1,m2,m3,k,:),[],slot) + reshape(Ivf(k,m1,ind(1),:),[],slot) + reshape(Ivf(k,m2,ind(2),:),[],slot);
                end
            end
            Ifv(k,m3,ind(3),:) = max(sIfv,[],1);
        end
        %Ifv(k,:,ind(3),:) = Ifv(k,:,ind(3),:) - Ifv(k,1,ind(3),:);
    end

    Ifv = Ifv - Ifv(:,1,:,:);

    
    
    %变量节点更新：Ivg
    for k =1:V*Nt
        ind = find(F(:,k)==1);
        Ifv_sum(1,:,k,:) = sum( Ifv(ind(1:len_col),:,k,:) );
        % analogue of normalization in MPA, it can be removed (s1 and s2), but at high SNR and/or number of iterations NaN LLR values can be exist, so Max-Log-MPA is required
        for i = 1:len_col
            Ivf(ind(i),:,k,:) = Ifv_sum(1,:,k,:) - Ifv(ind(i),:,k,:);
        end
    end

    %Ivf = log( exp(Ivf)./sum( exp(Ivf),2 ) );

end
%最后的符号概率：Q
Ifv_sum = Ifv_sum - max(Ifv_sum, [], 2);
P_Px = exp(permute(squeeze(Ifv_sum),[1 3 2]));
u_llr0 = Px2bitLLR(P_Px);
% u_llr0 = zeros(bit_num,slot,V*Nt);
% u_llr0(1,:,:) = log((P_Px(1,:,:) + P_Px(2,:,:))./(P_Px(3,:,:) + P_Px(4,:,:)));
% u_llr0(2,:,:) = log((P_Px(1,:,:) + P_Px(3,:,:))./(P_Px(2,:,:) + P_Px(4,:,:)));
u_llr = reshape(reshape(min( max(u_llr0,-10), 10),bit_num,slot*V*Nt),bit_num*slot*Nt,V)';

u_hat = (1 - sign(u_llr))/2;


end