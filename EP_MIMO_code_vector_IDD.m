function [u_llr,u_hat] = EP_MIMO_code_vector_IDD(y,K,V,M,F,CB,N0,h,Out_iter,Det_iter,Dec_iter,Nt,Nr,slot,bit_num,cfgLDPCDec,bit_table_normal,isinterleaved)

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
log_Px = zeros(M, V*Nt, slot);

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

Idec_R = zeros(bit_num*slot*Nt, V);
for k = 1 : Out_iter
    for kk = 1: Det_iter
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

        %I_F_to_V = exp(llr_I_F_to_V);

        for i = 1:V*Nt
            log_Px(:,i,:) = reshape(sum(llr_I_F_to_V(index_V(:,i),:,i,:)),M,1,slot);
            %llr(:,i,:) = squeeze(sum(llr_I_F_to_V(index_V(:,i),:,i,:)))';
        end
        log_Px = log_Px - max(log_Px, [], 1);
        Px = exp(log_Px);
        Px = Px ./ sum(Px);

        if(any(isnan(Px), 'all'))
            fprintf("det error!\n");
        end

        Mu = squeeze(sum(repmat(permute(Px,[4 1 2 3]),K*Nr,1,1,1).*CB,2));
        Sigma_temp = squeeze(sum(repmat(permute(Px,[4 1 2 3]),K*Nr,1,1,1).*abs(CB).^2,2)) - ...
            abs(Mu).^2;
        Sigma = max(Sigma_temp, 5e-7);

        %VN update

        Sigma_V_to_F_temp = (Sigma.*Sigma_F_to_V)./(Sigma_F_to_V - Sigma);
        Mu_V_to_F_temp = Sigma_V_to_F_temp.*(Mu./Sigma - Mu_F_to_V./Sigma_F_to_V);

        Sigma_V_to_F_index = Sigma_V_to_F_temp < 0;
        Sigma_V_to_F_temp(Sigma_V_to_F_index) = Sigma_V_to_F(Sigma_V_to_F_index);
        Mu_V_to_F_temp(Sigma_V_to_F_index) = Mu_V_to_F(Sigma_V_to_F_index);
        Sigma_V_to_F = Sigma_V_to_F_temp;
        Mu_V_to_F = Mu_V_to_F_temp;
    end
    scale = 1;
    P_Px = permute(Px,[1 3 2]);
    u_llr0 = Px2bitLLR(P_Px);
    Idet_P = reshape(u_llr0, bit_num*slot*Nt, V);
    Idet_R = scale * Idet_P;
    Idet_E = Idet_R;  

    if isinterleaved
        Idet_E_deinterleaved = reshape(permute(reshape(Idet_E, bit_num, size(Idet_E, 1) / bit_num, V), [2 1 3]), size(Idet_E, 1), V);
    else
        Idet_E_deinterleaved = Idet_E;
    end
    

    Idec_A = Idet_E_deinterleaved + Idec_R;
    Idec_P = ldpcDecode(Idec_A,cfgLDPCDec,Dec_iter,"DecisionType","soft","Termination","max","OutputFormat","whole");
    Idec_R = Idec_P - Idec_A;
    %Idec_R = min(max(Idec_R, -50), 50);
    mean_Idec_R = mean(abs(Idec_R), 1);
    mean_Idec_A = mean(abs(Idec_A), 1);
    Idec_E = Idec_R + 0 * mean_Idec_R./mean_Idec_A.*Idec_A;

    if isinterleaved
        Idct_E_interleaved = reshape(permute(reshape(Idec_E, size(Idec_E, 1) / bit_num, bit_num, V), [2 1 3]), size(Idec_E, 1), V);
    else
        Idct_E_interleaved = Idec_E;
    end

    Idet_A = Idct_E_interleaved ;

    Idec_E_reshape = reshape(Idet_A, bit_num, slot, Nt*V);

    llr_Idec_E_symbol = squeeze(sum(permute(-Idec_E_reshape,[1 4 2 3]).* bit_table_normal', 1));
    llr_Idec_E_symbol = llr_Idec_E_symbol - max(llr_Idec_E_symbol, [], 1);
    Px_Idec_E_symbol = permute(exp(llr_Idec_E_symbol) ./ sum(exp(llr_Idec_E_symbol), 1), [1 3 2]);

    Px = Px.*Px_Idec_E_symbol ./(sum(Px.*Px_Idec_E_symbol, 1));

    if(any(isnan(Px), 'all')) 
        fprintf("dec error!\n");
    end

    Mu = squeeze(sum(repmat(permute(Px,[4 1 2 3]),K*Nr,1,1,1).*CB,2));
    Sigma_temp = squeeze(sum(repmat(permute(Px,[4 1 2 3]),K*Nr,1,1,1).*abs(CB).^2,2)) - ...
        abs(Mu).^2;
    Sigma = max(Sigma_temp, 5e-7);
    Mu_V_to_F = Mu;
    Sigma_V_to_F = Sigma; 

end

u_llr = Idec_P;
u_hat = (1 - sign(u_llr))/2;

end