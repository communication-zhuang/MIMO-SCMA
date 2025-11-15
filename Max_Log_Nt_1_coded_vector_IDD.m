function [u_llr,u_hat] = Max_Log_Nt_1_coded_vector_IDD(y,K,V,M,F,CB,N0,h,Out_iter,Det_iter,Dec_iter,Nt,Nr,slot,bit_num,cfgLDPCDec,bit_table_normal, isinterleaved)
len_col = sum(F(:,1));
f = zeros(M, M, M, K*Nr,slot);
Ivf = zeros(K*Nr, M, V*Nt,slot);
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

%Idec_R = zeros(bit_num*slot*Nt, V);
for Out = 1 : Out_iter
    for Det = 1 : Det_iter
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
                Ifv(k,m1,ind(1),:) = max(sIfv, [], 1);
            end
            for m2 = 1:M
                sIfv = zeros(M*M, 1,slot);
                for m1 = 1:M
                    for m3 = 1:M
                        sIfv((m1-1)*M+m3,1,:) = f(m1,m2,m3,k,:) + Ivf0(k,m1,ind(1),1,:) + Ivf0(k,m3,ind(3),1,:);
                        %sIfv((m1-1)*M+m3,1,:) = reshape(f(m1,m2,m3,k,:),[],slot) + reshape(Ivf(k,m1,ind(1),:),[],slot) + reshape(Ivf(k,m3,ind(3),:),[],slot);
                    end
                end
                Ifv(k,m2,ind(2),:) = max(sIfv, [], 1);
            end
            for m3 = 1:M
                sIfv = zeros(M*M, 1,slot);
                for m1 = 1:M
                    for m2 = 1:M
                        sIfv((m1-1)*M+m2,1,:) = f(m1,m2,m3,k,:) + Ivf0(k,m1,ind(1),1,:) + Ivf0(k,m2,ind(2),1,:);
                        %sIfv((m1-1)*M+m2,1,:) = reshape(f(m1,m2,m3,k,:),[],slot) + reshape(Ivf(k,m1,ind(1),:),[],slot) + reshape(Ivf(k,m2,ind(2),:),[],slot);
                    end
                end
                Ifv(k,m3,ind(3),:) = max(sIfv, [], 1);
            end
        end

        Ifv = Ifv - Ifv(:,1,:,:);

        %变量节点更新：Ivg
        for k =1:V*Nt
            ind = find(F(:,k)==1);
            Ifv_sum(1,:,k,:) = sum( Ifv(ind(1:len_col),:,k,:) );
            % analogue of normalization in MPA, it can be removed (s1 and s2), but at high SNR and/or number of iterations NaN LLR values can be exist, so Max-Log-MPA is required
            for i = 1:len_col
                Ivf(ind(i),:,k,:) = Ifv_sum(1,:,k,:) - Ifv(ind(i),:,k,:);
                %Ivf(ind(i),:,k,:) = log( exp(Ivf(ind(i),:,k,:)) ./ sum(exp(Ivf(ind(i),:,k,:)), 2) );
            end
        end
    end
    
    Ifv_sum = Ifv_sum - max(Ifv_sum, [], 2);
    u = mean(Ifv_sum, [2 3]);
    MAD = mean(abs(Ifv_sum - u), [2 3]);
    gamma = 0.1;
    Ifv_sum = tanh(gamma * MAD) .* Ifv_sum;
    P_Px = exp(permute(squeeze(Ifv_sum),[1 3 2]));
    
    u_llr0 = Px2bitLLR(P_Px);
    Idet_P = reshape(u_llr0, bit_num*slot*Nt, V);
    Idet_R = Idet_P;
    Idet_E = Idet_R;

    if isinterleaved
        Idet_E_deinterleaved = reshape(permute(reshape(Idet_E, bit_num, size(Idet_E, 1) / bit_num, V), [2 1 3]), size(Idet_E, 1), V);
    else
        Idet_E_deinterleaved = Idet_E;
    end

    Idec_A = Idet_E_deinterleaved;
    Idec_P = ldpcDecode(Idec_A,cfgLDPCDec,Dec_iter,"DecisionType","soft","Termination","max","OutputFormat","whole");%"Termination","max",
    Idec_R = Idec_P - Idec_A;
    MI = mean(1 - log2 (exp(-abs(Idec_R)) + 1), 1);
    alpha = 10;
    theta = 0.4;
    a = 1 ./ (1 + exp(-alpha * (MI - theta)));
    %Idec_R = min(max(Idec_R, -50), 50);
    Idec_E = a .* Idec_R + (1 - a) .* Idec_A;
    %Idec_E = Idec_R + 0.6* mean_Idec_R./mean_Idec_A .* Idec_A;

    if isinterleaved
        Idct_E_interleaved = reshape(permute(reshape(Idec_E, size(Idec_E, 1) / bit_num, bit_num, V), [2 1 3]), size(Idec_E, 1), V);
    else
        Idct_E_interleaved = Idec_E;
    end


    Idet_A = Idct_E_interleaved;

    Idec_E_reshape = reshape(Idet_A, bit_num, slot, Nt*V);
    llr_Idec_E_symbol = squeeze(sum(permute(-Idec_E_reshape,[1 4 2 3]).* bit_table_normal', 1));
    Ivf = repmat(permute(llr_Idec_E_symbol, [4 1 3 2]), K*Nr, 1, 1, 1);
    %Ivf = Ivf1 ;%- Ifv;
    
end
u_llr = Idec_P;
u_hat = (1 - sign(u_llr))/2;
end