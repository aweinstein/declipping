% Compute the minimum requiered number of measurement requiered to recover a 
% signal, for different levels of sparsity.
% Figure 4 of the paper

clear all
close all
clc

N = 2^7;
K = 1:1:10;
methods = {
    'consOMP',
    'consOMP_realDFT',
    'tpcc'
};
cvx_quiet(true);
Ms = zeros(length(K), length(methods));
Cls = zeros(length(K), length(methods));
maxs = zeros(length(K),1);
cf = zeros(length(K), 1);
n_trials = 100;
tic
% For each sparsity level, we do a binary search to find the number of 
% measurements necesary for recovery.
k = 0;
for i = 1:length(K),
    for j = 1:length(methods),
        M_trial = zeros(n_trials, 1);
        cl_trial = zeros(n_trials, 1);
        for t = 1:n_trials,
            x = get_sparse_dct_signal(N, K(i));
            M = round(N / 2); % Initial guess
            high_M = N;
            low_M = 2;
            last_M = NaN;
            while 1,
                cl = get_clip_level(x, M);
                x_hat = declip_dct(x, cl, methods{j});
                err = norm(x - x_hat);
                if err > 1e-2, % Failed recovery, increase M 
                    low_M = M;
                    M = round((M + high_M) / 2);
                else           % Succesful recover, decrease M 
                    high_M = M;
                    M = round((M + low_M) / 2);
                end

                if last_M == M,
                    break
                end
                last_M = M;
            end
            M_trial(t) = M;
            cl_trial(t) = cl;
            
        end
        Ms(i,j) = mean(M_trial);
        Cls(i,j) = mean(cl_trial);
        k = k + 1;
        fprintf('Finish iteration %d out of %d at %s, for %s, K=%d \n', ...
                k, length(K)*length(methods), datestr(now,13), ...
                methods{j}, K(i))
    end
end
elapsed_time = toc

%% Save the results
save_results = true;
if save_results,
    file_name = [mfilename '_' datestr(now,30)]
    save(file_name, '-V7')
end



%% Plot the results
figure(1), clf
hold on
symb = 'oxd^sv';
labels = methods;
g = linspace(0,0.8,length(methods));
for j = 1:length(methods),
    p(j) = plot(K, Ms(:,j), [symb(j) '-']);
    set(p(j), 'LineWidth', 2, 'Color', g(j)*ones(3,1), 'MarkerSize', 7);
end
legend(methods)
% figure(1), clf
% % subplot(211) 
% hold on
% symb = 'oxd^sv';
% labels = {'ell1', 'rew ell1', 'OMP', 'OMP norm', 'LS-OMP', 'rand'};
% labels = methods;
% g = linspace(0,0.8,length(methods));
% for j = 1:length(methods),
%     p(j) = plot(K, Ms(:,j), [symb(j) '-']);
%     set(p(j), 'LineWidth', 2, 'Color', g(j)*ones(3,1), 'MarkerSize', 7);
% end
% 
% for j = 1:length(methods),
%     annotate(p(j), 0.25 + j, labels{j}, 'lr', pi/3, 0.03, 1, 14)
% end
% 
% xlabel('Sparsity K')
% ylabel('Mmin')
% 
% % subplot(212)
% % plot(K, cf,'bo-','linewidth',2)
% % xlabel('Sparsity K')
% % ylabel('Crest factor')
% % %legend('ell_1', 'rew ell_1', 'OMP', 'OMP norm', 'LS-OMP')
% savefig(['exp1_' datestr(now,30) '.pdf'],'pdf')
