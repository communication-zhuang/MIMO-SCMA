
clc;
clear;

rng(100);

tic
EbN0 = 5;

[CB,F,K,M,V] = get_default_CB(1); % K number of orthogonal resources; M number of codewords in each codebook; V number of users (layers)
iter_mpa_num = 8;

Nerrbits = zeros(1,length(EbN0));
NerrFrame = zeros(1,length(EbN0));

Nbits = zeros(1,length(EbN0));
BER   = zeros(1, length(EbN0));
FER = zeros(1,length(EbN0));


%maxNumErrs = 10000;
% maxNumBits = 1e7; %total numer of bits
% minNumBits = 50000;
minNumErrs = 100;

rate = 1;
bit_power = 1/6;
bit_num = log2(M);
Code=[0,0;0,1;1,0;1,1];

for iter_ebn0 = 1:length(EbN0)

    loop = 0;
    error_frame = 0;

    while (error_frame < minNumErrs) 
        loop = loop + 1;
        infobits = randi([0 1],V,bit_num);

        if (mod(loop,10) == 0)
            fprintf("\nNow Iter: %d\tNow SNR: %d\tNow Sigma: %f\tNow Error Frame: %d\tNow Error Bits: %d", ...
                loop, EbN0(iter_ebn0), N0, error_frame, Nerrbits(iter_ebn0));
        end
  
        x = bi2de(infobits,log2(M),'left-msb');
        h = sqrt(1/2)*(randn(K, V) + 1i*randn(K,V)); %Rayleigh channel
        %h = ones(K,V); % AWGN channel
        s = squeeze(scmaenc(x, CB, h)); 

        N0 = bit_power * K /(bit_num * rate * (10^(EbN0(iter_ebn0)/10)));
        sigma_sq = N0/2;
        noise = sqrt(sigma_sq)*(randn(size(s))+1i*randn(size(s))); 
        y = s + noise;
        
        %Factor graph calculation
        %u_hat = MaxLog_decode(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code);
        u_hat = LLR_Max_Log_decode(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code);
        %u_hat = Min_Max(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code);
        %u_hat = GAI_decode(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code);
        %u_hat = EP(y,K,V,M,F,CB,N0,h,iter_mpa_num,Code);
        
        
        %**********************************************************
        m_hat = reshape(u_hat',1,bit_num*V);
        m_reshape = reshape(infobits', 1, bit_num*V);  
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

