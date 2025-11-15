function u_hat = GAI_decode(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr)

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
Igv = zeros(K, V, M);
Pvg = 0.25*ones(K,V,M);
%LLR_Pri = zeros(M, V);
%P_Pri = exp(LLR_Pri - max(LLR_Pri));
%P_Pri = exp(LLR_Pri)./sum(exp(LLR_Pri));
%Pvg = repmat(permute(P_Pri,[3 2 1]),K,1);
CB0 = permute(CB,[1,3,2]);

%迭代开始

for mpa_iter = 1:iter_mpa_num

    %资源节点更新：Igv
    Mu_all = sum(Pvg.*CB0,3);
    Sigma_all = sum(Pvg.*abs(CB0).^2,3) - abs(sum(Pvg.*CB0,3)).^2;
    for k = 1:K
        ind = find(F(k,:)==1);
        for m1 = 1:M
            %             Igv(k,ind(1),m1) = -(1/N0)*abs( y(k) - (CB0(k,ind(1),m1)*h(k,ind(1))+Mu_all(k,ind(2))*h(k,ind(2))+Mu_all(k,ind(3))*h(k,ind(3)) ))^2 ...
            %                                +(1/N0)*abs( y(k) - (CB0(k,ind(1),1)*h(k,ind(1))+Mu_all(k,ind(2))*h(k,ind(2))+Mu_all(k,ind(3))*h(k,ind(3)) ))^2;
            Igv(k,ind(1),m1) = -(1/(N0+Sigma_all(k,ind(2))*abs(h(k,ind(2))).^2 + Sigma_all(k,ind(3))*abs(h(k,ind(3))).^2))*abs( y(k) - (CB0(k,ind(1),m1)*h(k,ind(1))+Mu_all(k,ind(2))*h(k,ind(2))+Mu_all(k,ind(3))*h(k,ind(3)) ))^2 ...
                +(1/(N0+Sigma_all(k,ind(2))*abs(h(k,ind(2))).^2 + Sigma_all(k,ind(3))*abs(h(k,ind(3))).^2))*abs( y(k) - (CB0(k,ind(1),1)*h(k,ind(1))+Mu_all(k,ind(2))*h(k,ind(2))+Mu_all(k,ind(3))*h(k,ind(3)) ))^2;
        end

        for m2 = 1:M
            Igv(k,ind(2),m2) = -(1/(N0+Sigma_all(k,ind(1))*abs(h(k,ind(1))).^2 + Sigma_all(k,ind(3))*abs(h(k,ind(3))).^2))*abs( y(k) - (Mu_all(k,ind(1))*h(k,ind(1))+CB0(k,ind(2),m2)*h(k,ind(2))+Mu_all(k,ind(3))*h(k,ind(3)) ))^2 ...
                +(1/(N0+Sigma_all(k,ind(1))*abs(h(k,ind(1))).^2 +  Sigma_all(k,ind(3))*abs(h(k,ind(3))).^2))*abs( y(k) - (Mu_all(k,ind(1))*h(k,ind(1))+CB0(k,ind(2),1)*h(k,ind(2))+Mu_all(k,ind(3))*h(k,ind(3)) ))^2;
        end

        for m3 = 1:M
            Igv(k,ind(3),m3) = -(1/(N0+Sigma_all(k,ind(1))*abs(h(k,ind(1))).^2 + Sigma_all(k,ind(2))*abs(h(k,ind(2))).^2))*abs( y(k) - (Mu_all(k,ind(1))*h(k,ind(1))+Mu_all(k,ind(2))*h(k,ind(2))+CB0(k,ind(3),m3)*h(k,ind(3)) ))^2 ...
                +(1/(N0+Sigma_all(k,ind(1))*abs(h(k,ind(1))).^2 + Sigma_all(k,ind(2))*abs(h(k,ind(2))).^2))*abs( y(k) - (Mu_all(k,ind(1))*h(k,ind(1))+Mu_all(k,ind(2))*h(k,ind(2))+CB0(k,ind(3),1)*h(k,ind(3)) ))^2;
        end

    end



    %变量节点更新：Ivg
    for k =1:V
        ind = find(F(:,k)==1);
        Ivg(ind(1),k,:) = Igv(ind(2),k,:);
        Ivg(ind(2),k,:) = Igv(ind(1),k,:);

    end
    Pvg = exp(Ivg)./sum(exp(Ivg),3);
    
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