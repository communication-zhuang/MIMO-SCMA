function u_hat = LLR_Max_Log_decode(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code)


f = zeros(M, M, M, K);

for k = 1:K % resourses
    ind = find(F(k,:)==1); % non-zero elements, paths
    for m1 = 1:M
        for m2 = 1:M
            for m3 = 1:M
                f(m1,m2,m3,k) = -(1/N0)*abs(y(k)-(CB(k,m1,ind(1))*h(k,ind(1))+CB(k,m2,ind(2))*h(k,ind(2))+CB(k,m3,ind(3))*h(k,ind(3))  ))^2;
            end
        end
    end
end

Ivg = zeros(K, V, M);
%Ivg = log(1/M*ones(K, V, M));
Igv = zeros(K, V, M);


%迭代开始
for mpa_iter = 1:iter_mpa_num

    %资源节点更新：Igv
    for k = 1:K
        ind = find(F(k,:)==1);
        for m1 = 1:M
            sIgv = zeros(M*M,1);
            for m2 = 1:M
                for m3 = 1:M
                    sIgv((m2-1)*M+m3,1) = f(m1,m2,m3,k) + Ivg(k,ind(2),m2) + Ivg(k,ind(3),m3);
                end
            end
            Igv(k,ind(1),m1) = max(sIgv);
        end

        for m2 = 1:M
            sIgv = zeros(M*M, 1);
            for m1 = 1:M
                for m3 = 1:M
                    sIgv((m1-1)*M+m3,1) = f(m1,m2,m3,k) + Ivg(k,ind(1),m1) + Ivg(k,ind(3),m3);
                end
            end
            Igv(k,ind(2),m2) = max(sIgv);
        end

        for m3 = 1:M
            sIgv = zeros(M*M, 1);
            for m1 = 1:M
                for m2 = 1:M
                    sIgv((m1-1)*M+m2,1) = f(m1,m2,m3,k,:) + Ivg(k,ind(1),m1) + Ivg(k,ind(2),m2);
                end
            end
            Igv(k,ind(3),m3) = max(sIgv);
        end
    end
    Igv(:, :, :) = Igv(:, :, :) - Igv(:, :, 1);



    %变量节点更新：Ivg
    for k =1:V
        ind = find(F(:,k)==1);
         Ivg(ind(1),k,:) = Igv(ind(2),k,:) ;
         Ivg(ind(2),k,:) = Igv(ind(1),k,:) ;
        
    end

   
end



%最后的符号概率：Q
Q = zeros(M, V);
for k = 1:V
    ind = find(F(:,k)==1);
    for m = 1:M
        Q(m,k) = Igv(ind(1),k,m) + Igv(ind(2),k,m);
    end
end

% Q = log(exp(Q)./sum(exp(Q)));
u_hat = zeros(V,log2(M));
for k=1:V
    [~,L] = max(Q(:,k));
    u_hat(k,:) = Code(L,:);
end
end