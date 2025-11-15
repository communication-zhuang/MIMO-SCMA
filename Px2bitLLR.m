function u_llr0 = Px2bitLLR(Px)

bit_num = log2(size(Px, 1));
u_llr0 = zeros(bit_num, prod(size(Px))/size(Px, 1));
thre_min = 1e-10;
iscut = 1;

switch bit_num
    case 2
        if iscut
            bit1_0to1 = min( max((Px(1,:) + Px(2,:)) ./ (Px(3,:) + Px(4,:)), thre_min / (1 - thre_min)), (1 - thre_min) / thre_min);
            bit2_0to1 = min( max((Px(1,:) + Px(3,:)) ./ (Px(2,:) + Px(4,:)), thre_min / (1 - thre_min)), (1 - thre_min) / thre_min);
            u_llr0(1,:) = log(bit1_0to1);
            u_llr0(2,:) = log(bit2_0to1);
        else
            u_llr0(1,:) = log((Px(1,:) + Px(2,:)) ./ (Px(3,:) + Px(4,:)));
            u_llr0(2,:) = log((Px(1,:) + Px(3,:)) ./ (Px(2,:) + Px(4,:)));
        end


    case 3
        bit1_0to1 = min( max(sum(Px(1:4, :)) ./ sum(Px(5:8, :)), thre_min / (1 - thre_min)), (1 - thre_min) / thre_min);
        bit2_0to1 = min( max(sum(Px([1:2 5:6], :)) ./ sum(Px([3:4 7:8], :)), thre_min / (1 - thre_min)), (1 - thre_min) / thre_min);
        bit3_0to1 = min( max(sum(Px([1 3 5 7], :)) ./ sum(Px([2 4 6 8], :)), thre_min / (1 - thre_min)), (1 - thre_min) / thre_min);
        u_llr0(1,:) = log(bit1_0to1);
        u_llr0(2,:) = log(bit2_0to1);
        u_llr0(3,:) = log(bit3_0to1);
end

end
