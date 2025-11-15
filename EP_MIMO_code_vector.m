function [u_llr,u_hat] = EP_MIMO_code_vector(y,K,V,M,F,CB,N0,h,iter_num,Nt,Nr,slot,bit_num)

Mu_V_to_F = zeros(K*Nr,V*Nt,slot);
Sigma_V_to_F = zeros(K*Nr,V*Nt,slot);
Mu_F_to_V = zeros(K*Nr,V*Nt,slot);
Sigma_F_to_V = zeros(K*Nr,V*Nt,slot);
I_F_to_V = zeros(K*Nr,M,V*Nt,slot);
Mu_all = zeros(K*Nr,1,slot);
Sigma_all = zeros(K*Nr,1,slot);
Mu = zeros(K*Nr,V*Nt,slot);
Sigma = zeros(K*Nr,V*Nt,slot);
Px = zeros(M,V*Nt,slot);

%preprocess
for i = 1:size(F,1)
    index_F(i,:) = find(F(i,:)==1);
end
for i = 1:size(F,2)
    index_V(:,i) = find(F(:,i)==1);
end

Max = 100;

%Initialization
Mu_V_to_F = zeros(K*Nr,V*Nt,slot);
Sigma_V_to_F = Max*ones(K*Nr,V*Nt,slot);
Px = 1/M*ones(M,V*Nt,slot);


for iter = 1:iter_num
    %FN update
    Mu_H = h.*Mu_V_to_F;
    Sigama_H = (abs(h).^2).*Sigma_V_to_F;

    for i = 1:K*Nr
        Mu_all(i,1,:) = sum(Mu_H(i,index_F(i,:),:),2);
        Sigma_all(i,1,:) = sum(Sigama_H(i,index_F(i,:),:),2);
    end
    Mu_F_to_V = (permute(y,[1 3 2]) - (Mu_all - Mu_H))./h;
    Sigma_F_to_V = (Sigma_all - Sigama_H + N0)./abs(h).^2;

    llr_I_F_to_V = -abs(CB - repmat(permute(Mu_F_to_V,[1 4 2 3]),1,M,1,1)).^2./ ...
              repmat(permute(Sigma_F_to_V,[1 4 2 3]),1,M,1,1);

    I_F_to_V = exp(llr_I_F_to_V);
    
    for i = 1:V*Nt
        Px(:,i,:) = reshape(prod(I_F_to_V(index_V(:,i),:,i,:)),M,1,slot);
        %llr(:,i,:) = squeeze(sum(llr_I_F_to_V(index_V(:,i),:,i,:)))';
    end
    Px = Px./sum(Px);

    if(any(isnan(Px), 'all'))
        fprintf("Px error!\n");
    end

    Mu = squeeze(sum(repmat(permute(Px,[4 1 2 3]),K*Nr,1,1,1).*CB,2));
    Sigma_temp = squeeze(sum(repmat(permute(Px,[4 1 2 3]),K*Nr,1,1,1).*abs(CB).^2,2)) - ...
        abs(Mu).^2;
    Sigma_index = Sigma_temp < 5e-7;
    Sigma_temp(Sigma_index) = 5e-7;
    Sigma = Sigma_temp;
    
    %VN update
    
    Sigma_V_to_F_temp = (Sigma.*Sigma_F_to_V)./(Sigma_F_to_V - Sigma);
    Mu_V_to_F_temp = Sigma_V_to_F_temp.*(Mu./Sigma - Mu_F_to_V./Sigma_F_to_V);
    
    Sigma_V_to_F_index = Sigma_V_to_F_temp < 0;
    Sigma_V_to_F_temp(Sigma_V_to_F_index) = Sigma_V_to_F(Sigma_V_to_F_index);
    Mu_V_to_F_temp(Sigma_V_to_F_index) = Mu_V_to_F(Sigma_V_to_F_index);
    Sigma_V_to_F = Sigma_V_to_F_temp;
    Mu_V_to_F = Mu_V_to_F_temp;
end

P_Px = permute(Px,[1 3 2]);
u_llr0 = Px2bitLLR(P_Px);
% u_llr0 = zeros(bit_num,slot,V*Nt);
% u_llr0(1,:,:) = log((P_Px(1,:,:) + P_Px(2,:,:))./(P_Px(3,:,:) + P_Px(4,:,:)));
% u_llr0(2,:,:) = log((P_Px(1,:,:) + P_Px(3,:,:))./(P_Px(2,:,:) + P_Px(4,:,:)));
u_llr = reshape(reshape(min( max(u_llr0,-10), 10),bit_num,slot*V*Nt),bit_num*slot*Nt,V)';
u_hat = (1 - sign(u_llr))/2;

end