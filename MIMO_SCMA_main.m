clc;
clear;

rng(100);
tic
EbN0 = -2; 
iter_mpa_num = 6;

CodewordType = 1; % 1. huawei  2. AMICB1 3. AMICB2

[CB0,F0,K,M,V] = get_default_CB(CodewordType); % K number of orthogonal resources; M number of codewords in each codebook; V number of users (layers)
Nt = 1;
Nr = 4;


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


Nerrbits = zeros(1,length(EbN0));
NerrFrame = zeros(1,length(EbN0));

Nbits = zeros(1,length(EbN0));
BER   = zeros(1, length(EbN0));
FER = zeros(1,length(EbN0));


%maxNumErrs = 10000;
% maxNumBits = 1e7; %total numer of bits
% minNumBits = 50000;
minNumErrs = 200;


Es = 1/6;
bit_num = log2(M);
RE_rate = bit_num / K;
Code_rate = 1;
switch M
    case 4 
        Code=[0,0;0,1;1,0;1,1];
    case 8
        Code = [0,0,0;0,0,1;0,1,0;0,1,1;1,0,0;1,0,1;1,1,0;1,1,1];
end


for iter_ebn0 = 1:length(EbN0)

    loop = 0;
    error_frame = 0;

    while (error_frame < minNumErrs) 
        loop = loop + 1;
        infobits = randi([0 1],V*Nt,bit_num);

        if (mod(loop,10) == 0)
            fprintf("\nNow Iter: %d\tNow SNR: %d\tNow Sigma: %f\tNow Error Frame: %d\tNow Error Bits: %d", ...
                loop, EbN0(iter_ebn0), N0, error_frame, Nerrbits(iter_ebn0));
        end
  
        x = bi2de(infobits,'left-msb');
        %h = sqrt(1/2)*(randn(K*Nr, V*Nt) + 1i*randn(K*Nr,V*Nt)); %Rayleigh channel
        h = ones(K*Nr,V*Nt); % AWGN channel
        s = squeeze(scmaenc(x, CB, h)); 

        N0 = Es /(RE_rate * Code_rate * (10^(EbN0(iter_ebn0)/10)));
        sigma_sq = N0/2;
        noise = sqrt(sigma_sq)*(randn(size(s))+1i*randn(size(s))); 
        y = s + noise;
        
        %Factor graph calculation
        %u_hat = MaxLog_decode_MIMO(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr,index_MPA);
        u_hat = BsMP_decode_Nt_1(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr);
        %[u_llr, u_hat] = EP_MIMO(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr);
        %u_hat = BsMP_MIMO(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code,Nt,Nr,index_MPA);
        
        %**********************************************************
        m_hat = reshape(u_hat',1,bit_num*V*Nt);
        m_reshape = reshape(infobits', 1, bit_num*V*Nt);  
        err = sum(m_hat~=m_reshape);
        if err ~= 0
            error_frame = error_frame + 1;
        end
        Nerrbits(iter_ebn0) = Nerrbits(iter_ebn0) + err;
        Nbits(iter_ebn0) = Nbits(iter_ebn0) + length(m_reshape);   
    end
    

    BER(iter_ebn0) = Nerrbits(iter_ebn0)/Nbits(iter_ebn0);	
    FER(iter_ebn0) = error_frame/loop;
	
        
end

toc

semilogy(EbN0,BER);

