clear all
close all
clc

N = 2^7;
K = [1:10];
methods = {
    'consOMP',
    'tpcc',
};
cvx_quiet(true);

M = 80;
n_trials = 1000;
opts.A = true;
opts.f = true; 
opts.phase = false;

n_tot = length(K) * n_trials * length(methods);
k = 1;
result = zeros(length(K), n_trials, length(methods), 5);
tic
for i = 1:length(K),
    for n = 1:n_trials,
        x = get_sparse_dct_signal(N, K(i));
        cl = get_clip_level(x, M);
        x_c = clip_signal(x, cl);   
        for m = 1:length(methods),
            x_hat = declip_dct(x, cl, methods{m});
            disp(['Iteration ' num2str(k) ' of ' num2str(n_tot)])
            k = k + 1;
            result(i,n,m,1) = norm(x - x_hat);
            result(i,n,m,2) = cl;
            result(i,n,m,3) = M;
        end
    end
end
elapsed_time = toc
%% Compute the probability of recovery (number of recoveries / n_trials
p_recovery = zeros(length(methods), length(K));
for m = 1:length(methods),
    for i = 1:length(K),
        err = result(i,:,m,1);
        p_recovery(m,i) = sum(err<1e-3) / n_trials;
    end
end

%% Save the results
save_results = false;
if save_results,
    file_name = [mfilename '_' datestr(now,30)];
    save(file_name, '-V7')
    datestr(now)
end

%% Plot the results
figure
plot(K, p_recovery(1,:), 'bx-')
hold on
plot(K, p_recovery(2,:), 'ro-')
legend('consOMP', 'TCPP')
ylim([0, 1.1])
