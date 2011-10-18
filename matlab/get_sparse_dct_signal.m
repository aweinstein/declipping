function x = get_sparse_dct_signal(N, K)

D = idct(eye(N));
rp = randperm(N);
rp(rp==1) = []; % We don't want to choose the DC term
support = rp(1:K);
A = sign(randn(K,1)) .* (0.5+rand(K,1));
x = D(:, support) * A;
