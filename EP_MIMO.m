function [u_llr,u_hat] = EP_MIMO(y,K,V,M,F,CB,N0,h,iter_num,Code,Nt,Nr)

Mu_V_to_F = zeros(K*Nr,V*Nt);
Sigma_V_to_F = zeros(K*Nr,V*Nt);
Mu_F_to_V = zeros(K*Nr,V*Nt);
Sigma_F_to_V = zeros(K*Nr,V*Nt);
I_F_to_V = zeros(K*Nr,M,V*Nt);
Mu_all = zeros(K*Nr,1);
Sigma_all = zeros(K*Nr,1);
Mu = zeros(K*Nr,V*Nt);
Sigma = zeros(K*Nr,V*Nt);
Px = zeros(M,V*Nt);

%preprocess
for i = 1:size(F,1)
    index_F(i,:) = find(F(i,:)==1);
end
for i = 1:size(F,2)
    index_V(:,i) = find(F(:,i)==1);
end

Max = 100;

%Initialization
Mu_V_to_F = zeros(K*Nr,V*Nt);
Sigma_V_to_F = Max*ones(K*Nr,V*Nt);
Px = 1/M*ones(M,V*Nt);


for iter = 1:iter_num
    %FN update
    Mu_H = h.*Mu_V_to_F;
    Sigama_H = (abs(h).^2).*Sigma_V_to_F;

    for i = 1:K*Nr
        Mu_all(i,1) = sum(Mu_H(i,index_F(i,:)));
        Sigma_all(i,1) = sum(Sigama_H(i,index_F(i,:)));
    end
    Mu_F_to_V = (y - (Mu_all - Mu_H))./h;
    Sigma_F_to_V = (Sigma_all - Sigama_H + N0)./abs(h).^2;

    llr_I_F_to_V = -abs(CB - repmat(permute(Mu_F_to_V,[1 3 2]),1,M,1)).^2./ ...
              repmat(permute(Sigma_F_to_V,[1 3 2]),1,M,1);
    I_F_to_V = exp(llr_I_F_to_V);

    for i = 1:V*Nt
        Px(:,i) = prod(I_F_to_V(index_V(:,i),:,i))';
       %llr(:,i) = sum(llr_I_F_to_V(index_V(:,i),:,i))';   
    end
    Px = Px./sum(Px);

    Mu = squeeze(sum(repmat(permute(Px,[3 1 2]),K*Nr,1,1).*CB,2));
    Sigma_temp = squeeze(sum(repmat(permute(Px,[3 1 2]),K*Nr,1,1).*abs(CB).^2,2)) - ...
        abs(Mu).^2;
    Sigma_index = Sigma_temp < 5e-17;
    Sigma_temp(Sigma_index) = 5e-17;
    Sigma = Sigma_temp;
    
    %VN update
%     Sigma_V_to_F = (Sigma.*Sigma_F_to_V)./(Sigma_F_to_V - Sigma);
%     Mu_V_to_F = Sigma_V_to_F.*(Mu./Sigma - Mu_F_to_V./Sigma_F_to_V);
    
    Sigma_V_to_F_temp = (Sigma.*Sigma_F_to_V)./(Sigma_F_to_V - Sigma);
    Mu_V_to_F_temp = Sigma_V_to_F_temp.*(Mu./Sigma - Mu_F_to_V./Sigma_F_to_V);
    
    Sigma_V_to_F_index = Sigma_V_to_F_temp < 0;
    Sigma_V_to_F_temp(Sigma_V_to_F_index) = Sigma_V_to_F(Sigma_V_to_F_index);
    Mu_V_to_F_temp(Sigma_V_to_F_index) = Mu_V_to_F(Sigma_V_to_F_index);
    Sigma_V_to_F = Sigma_V_to_F_temp;
    Mu_V_to_F = Mu_V_to_F_temp;
end

%u_llr0 = [( max(llr(1:2,:)) - max(llr(3:4,:)) )' ( max(llr([1 3],:)) - max(llr([2 4],:)) )'];
u_llr0 = Px2bitLLR(Px)';
%u_llr0 = [ log((Px(1,:) + Px(2,:))./(Px(3,:) + Px(4,:))) ; log((Px(1,:) + Px(3,:))./(Px(2,:) + Px(4,:)))  ]';
u_llr = min( max(u_llr0,-3), 3);

u_hat = (1 - sign(u_llr))/2;

% u_hat = zeros(V*Nt,log2(M));
% 
% for k=1:V*Nt
%     [~,L] = max(Px(:,k));
%     u_hat(k,:) = Code(L,:);
% end

end