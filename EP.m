function u_hat = EP(y,K,V,M,F,CB,N0,h,iter_num,Code)

Mu_V_to_F = zeros(K,V);
Sigma_V_to_F = zeros(K,V);
Mu_F_to_V = zeros(K,V);
Sigma_F_to_V = zeros(K,V);
I_F_to_V = zeros(K,M,V);
Mu_all = zeros(K,1);
Sigma_all = zeros(K,1);
Mu = zeros(K,V);
Sigma = zeros(K,V);
Px = zeros(M,V);

%preprocess
for i = 1:size(F,1)
    index_F(i,:) = find(F(i,:)==1);
end
for i = 1:size(F,2)
    index_V(:,i) = find(F(:,i)==1);
end

Max = 10000;

%Initialization
Mu_V_to_F = zeros(K,V);
Sigma_V_to_F = Max*ones(K,V);
Px = 1/M*ones(M,V);


for iter = 1:iter_num
    %FN update
    Mu_H = h.*Mu_V_to_F;
    Sigama_H = (abs(h).^2).*Sigma_V_to_F;

    for i = 1:K
        Mu_all(i,1) = sum(Mu_H(i,index_F(i,:)));
        Sigma_all(i,1) = sum(Sigama_H(i,index_F(i,:)));
    end
    Mu_F_to_V = (y - (Mu_all - Mu_H))./h;
    Sigma_F_to_V = (Sigma_all - Sigama_H + N0)./abs(h).^2;

%     for i = 1:V
%         if sum(Sigma_F_to_V(index_V(:,i),i)) <0.1
%             Sigma_F_to_V(:,i) = 0.1*Sigma_F_to_V(:,i)./sum(Sigma_F_to_V(index_V(:,i),i));
%         end
%     end

    I_F_to_V = exp(-abs(CB - repmat(permute(Mu_F_to_V,[1 3 2]),1,M,1)).^2./ ...
              repmat(permute(Sigma_F_to_V,[1 3 2]),1,M,1));
    

    for i = 1:V
        Px(:,i) = prod(I_F_to_V(index_V(:,i),:,i))';
    end
    Px = Px./sum(Px);

    Mu = squeeze(sum(repmat(permute(Px,[3 1 2]),K,1,1).*CB,2));
    Sigma_temp = squeeze(sum(repmat(permute(Px,[3 1 2]),K,1,1).*abs(CB).^2,2)) - ...
        abs(Mu).^2;
    Sigma_index = Sigma_temp < 5e-7;
    Sigma_temp(Sigma_index) = 5e-7;
    Sigma = Sigma_temp;
    
    %VN update
    Sigma_V_to_F = (Sigma.*Sigma_F_to_V)./(Sigma_F_to_V - Sigma);
    Mu_V_to_F = Sigma_V_to_F.*(Mu./Sigma - Mu_F_to_V./Sigma_F_to_V);
end

u_hat = zeros(V,log2(M));
for k=1:V
    [~,L] = max(Px(:,k));
    u_hat(k,:) = Code(L,:);
end





end