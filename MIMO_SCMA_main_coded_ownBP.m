
clc;
clear;

tic
EbN0 = 5; 
iter_mpa_num=3;
minNumErrs = 10;

[CB0,F0,K,M,V] = get_default_CB(); % K number of orthogonal resources; M number of codewords in each codebook; V number of users (layers)
Nt = 1;
Nr = 2;
RE_rate = 1/2;
Es = 1/6;
bit_num = log2(M);
Code=[0,0;0,1;1,0;1,1];

%LDPC Parameter only support code_N = 648 1294 1944;code_R = 1/2 2/3 3/4 5/6
code_N = 648;
code_R = 3/4;
code_K = floor(code_N * code_R);
code_M = code_N - code_K;
LDPCIter = 20;

%Get H and G  
if(exist(['matrix/H',num2str(code_M),'x',num2str(code_N),'.mat']) && exist(['matrix/G',num2str(code_K),'x',num2str(code_N),'.mat']))
    load(['matrix/H',num2str(code_M),'x',num2str(code_N)],'H');
    load(['matrix/G',num2str(code_K),'x',num2str(code_N)],'G');
else
    H = GetLDPC_H(code_N,code_R);
    G = GetLDPC_G(H);
    save(['matrix/H',num2str(code_M),'x',num2str(code_N)],'H');
    save(['matrix/G',num2str(code_K),'x',num2str(code_N)],'G');
end






for n1 = 1:K  %Extend Factor Graph
    F_row((n1-1)*Nr+1:n1*Nr,:) = repmat(F0(n1,:),Nr,1);
end
for n2 = 1:V
    F_col(:,(n2-1)*Nt+1:n2*Nt) = repmat(F_row(:,n2),1,Nt);
end
F = F_col;  


for n1 = 1:K   %Extend Codebok CB
    CB_row((n1-1)*Nr+1:n1*Nr,:,:) = repmat(CB0(n1,:,:),Nr,1,1);
end
for n2 = 1:V
    CB_col(:,:,(n2-1)*Nt+1:n2*Nt) = repmat(CB_row(:,:,n2),1,1,Nt);
end
CB = CB_col;  


% Preprocessing
len_row = sum(F(1,:));
len_col = sum(F(:,1));

index_M = zeros(M^len_row,len_row);
for i = 1:M^len_row
    for j = 1:len_row
        index_M(i,j) = mod(floor( (i-1)/(M^(j-1)) ), M) + 1;
    end
end

index_MPA = zeros(K*Nr,len_row,M^len_row);
for k = 1:K*Nr            
    ind = find(F(k,:)==1);
    for i = 1:M^(len_row)
        index_MPA(k,:,i) = sub2ind(size(CB),k*ones(1,len_row),index_M(i,:),ind(1:len_row));
    end
end


Nerrbits_uncoded = zeros(1,length(EbN0));
Nerrbits_coded = zeros(1,length(EbN0));
NerrFrame = zeros(1,length(EbN0));

Nbits = zeros(1,length(EbN0));
BER   = zeros(1, length(EbN0));
FER = zeros(1,length(EbN0));




%maxNumErrs = 10000;
% maxNumBits = 1e7; %total numer of bits
% minNumBits = 50000;


slot = code_N/(bit_num*Nt);
u_hat = zeros(V*Nt,slot*bit_num);
u_llr = zeros(V*Nt,slot*bit_num);



for iter_ebn0 = 1:length(EbN0)

    loop = 0;
    error_frame = 0;

    while (error_frame < minNumErrs) 
        loop = loop + 1;
        uncoded_bits = randi([0 1],V, code_K);
%         coded_bits = reshape([uncoded_bits mod(uncoded_bits*G(:,code_K+1:code_N),2)]', bit_num, slot*Nt*V)';
%         sym_index = reshape(bi2de(coded_bits,'left-msb'),[slot,Nt*V])';
        coded_bits = [uncoded_bits mod(uncoded_bits*G(:,code_K+1:code_N),2)];
        sym_index = reshape(bit2int(coded_bits',bit_num),[slot,Nt*V])';
        
        
        %P_coded_bits = reshape(permute(coded_bits,[2 1 3]),[V*Nt*bit_num code_N/(bit_num*Nt)]);
   
        if (mod(loop,10) == 0)
            fprintf("\nNow Iter: %d\tNow SNR: %d\tNow Sigma: %f\tNow Error Frame: %d\tNow Error Bits: %d", ...
                loop, EbN0(iter_ebn0), N0, error_frame, Nerrbits_coded(iter_ebn0));
        end
  
        %x = bi2de(infobits,log2(M),'left-msb');
        h = sqrt(1/2)*(randn(K*Nr, V*Nt) + 1i*randn(K*Nr,V*Nt)); %Rayleigh channel
        %h = ones(K*Nr,V*Nt); % AWGN channel
        s = scmaenc(sym_index, CB, h); 

        N0 = Es/(RE_rate*(10^(EbN0(iter_ebn0)/10)));
        sigma_sq = N0/2;
        noise = sqrt(sigma_sq)*(randn(size(s))+1i*randn(size(s))); 
        y = s + noise;
        
        %Factor graph calculation
        %u_hat = MaxLog_decode_MIMO(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr,index_MPA);
        %u_hat = BsMP_decode_Nt_1(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr);
        
        for i = 1:slot
            [u_llr0(:,(i-1)*bit_num+1:i*bit_num),u_hat0(:,(i-1)*bit_num+1:i*bit_num)] = EP_MIMO(y(:,i),K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr);
        end

        u_llr = reshape(u_llr0',code_N,V);

        [c_hat,c_llr] = bp_LDPC_decode(u_llr,H,V,code_N,code_M,LDPCIter);

        


        
        
        %**********************************************************
        m_hat = reshape(u_hat0',1,[]);
        k_hat = reshape(c_hat(1:code_K,:),1,[]);

        m_reshape = reshape(coded_bits', 1, []); 
        k_reshape = reshape(uncoded_bits',1,[]);
        err_uncoded = sum(m_hat ~= m_reshape);
        err_coded = sum(k_hat ~= k_reshape);
        if err_coded ~= 0
            error_frame = error_frame + 1;
        end

        Nerrbits_uncoded(iter_ebn0) = Nerrbits_uncoded(iter_ebn0) + err_uncoded;
        Nerrbits_coded(iter_ebn0) = Nerrbits_coded(iter_ebn0) + err_coded;
        
    end
    Nbits_uncoded(iter_ebn0) = loop * V * code_N; 
    Nbits_coded(iter_ebn0) = loop * V * code_K;  
    

    BER_uncoded(iter_ebn0) = Nerrbits_uncoded(iter_ebn0)/Nbits_uncoded(iter_ebn0);	
    BER_coded(iter_ebn0) = Nerrbits_coded(iter_ebn0)/Nbits_coded(iter_ebn0);	
    %FER(iter_ebn0) = error_frame/loop;
	
        
end

toc

semilogy(EbN0,BER);

