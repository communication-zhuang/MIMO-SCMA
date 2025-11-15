function G = GetLDPC_G(H)

[m,n] = size(H);
k = n - m;

gfH = gf(H);
Hs = gfH(:,1:k);
Hp = gfH(:,k+1:n);
Hp_inv = inv(Hp);            % Very slow operation
G = [eye(k) (Hp_inv*Hs).'];
G = double(G.x);


end