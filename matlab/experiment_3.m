% For a fixed method, compute probability of recovery as a function of sparsity
% level. Repeat for different values of number of measurements

clear all
close all
clc

N = 2^7;
K = [2:2:20];    % Sparsity level
method = 'tpcc';
cvx_quiet(true);

M = 50:20:90; % Number of measurements
n_trials = 1000;
opts.A = true;
opts.f = true; 
opts.phase = true;

n_tot = length(K) * n_trials * length(M);
k = 1;
result = zeros(length(K), n_trials, length(M), 5);
tic
for m = 1:length(M),
    for i = 1:length(K),
        for n = 1:n_trials,
            x = get_sparse_signal(N,1,K(i),[],[],[],opts);
            cl = get_clip_level(x, M(m));
            x_c = clip_signal(x, cl);   
            x_hat = declip(x, cl, method);
            if rem(k,100) == 0,
                disp(['Iteration ' num2str(k) ' of ' num2str(n_tot)])
            end
            k = k + 1;
            result(i,n,m,1) = norm(x - x_hat);
            result(i,n,m,2) = cl;
            result(i,n,m,3) = M(m);
            %result(i,n,m,4) = n_bumbs(x);
            %result(i,n,m,5) = n_bumbs(x_c);
        end
    end
end
elapsed_time = toc

%% Compute the probability of recovery (number of recoveries / n_trials
p_recovery = zeros(length(M), length(K));
for m = 1:length(M),
    for i = 1:length(K),
        err = result(i,:,m,1);
        p_recovery(m,i) = sum(err<1e-3) / n_trials;
    end
end

%% Save the results
save_results = true;
if save_results,
    file_name = [mfilename '_' datestr(now,30)];
    save(file_name, '-V7')
    datestr(now)
end

plot_results = false;
if plot_results,
    %% Plot the results
    figure(1), clf
    for i = 1:length(M),
        plot(K, p_recovery(i,:), 'bx-')
        hold on
    end
    ylim([0, 1.1])
end
