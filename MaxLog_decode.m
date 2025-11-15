function u_hat = MaxLog_decode(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code)

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

Ivg = log(1/M*ones(K, V, M));
Igv = zeros(K, V, M);

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


    %变量节点更新：Ivg
    for k =1:V
        ind = find(F(:,k)==1);
        %             s1 = log(sum(exp(Igv(ind1(1),k,:,:)),3));
        %             s2 = log(sum(exp(Igv(ind1(2),k,:,:)),3));

        % analogue of normalization in MPA, it can be removed (s1 and s2), but at high SNR and/or number of iterations NaN LLR values can be exist, so Max-Log-MPA is required
        Ivg(ind(1),k,:) = Igv(ind(2),k,:); %- repmat(s2,1,1,M,1);
        Ivg(ind(1),k,:) = log( exp( Ivg( ind(1),k,: ) )./sum( exp( Ivg( ind(1),k,: ) ) ) );

        Ivg(ind(2),k,:) = Igv(ind(1),k,:);% - repmat(s1,1,1,M,1);
        Ivg(ind(2),k,:) = log(exp(Ivg(ind(2),k,:))./sum(exp(Ivg(ind(2),k,:))));
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