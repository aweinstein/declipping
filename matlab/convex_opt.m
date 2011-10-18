% Decliping a sinusoid using a BP, BPCC and Reweighted ell_1 with CC
% Figure 1 of the paper

close all
clear all
clc

cvx_solver SDPT3 % CVX fails when using SeDuMi
cvx_quiet true;
%%
N = 2^7;        
n = (0:N-1)';
T = N;
K = 1;
M = 70; % Recovery
M = 66; % BP and BPCC fail
x = get_sparse_signal(N,1,K,[1],[1],[pi/4]);

cl = get_clip_level(x, M);
i_uc = x >= cl;
i_lc = x <= -cl;
x_c = x;
x_c(i_uc) = cl;
x_c(i_lc) = -cl;

%% Take as measurements the values of x that are not clipped
I = eye(N);
Phi = I(~(i_uc | i_lc),:);
y = Phi * x_c;

% Sparsity basis
Psi = ifft(eye(N));

% Solve using CVX, without the clipping constraints
disp('Solving BP ...')
cvx_begin
    variable alpha_hat(N) complex;
    minimize(norm(alpha_hat,1));
    subject to
        Phi * Psi * alpha_hat == y
cvx_end
x_hat_bp = real(ifft(alpha_hat));
%alpha_hat_l1_wo = alpha_hat;

% Solve using CVX, with the clipping constraints
disp('Solving BPCC ...')
Psi_u = Psi(i_uc,:);
Psi_l = Psi(i_lc,:);
cvx_begin
    variable alpha_hat(N) complex;
    minimize(norm(alpha_hat,1));
    subject to
        Phi * Psi * alpha_hat == y;
        real(Psi_u * alpha_hat) >= cl;
        real(Psi_l * alpha_hat) <= -cl;
cvx_end

x_hat_bpcc = real(ifft(alpha_hat));
alpha_hat_l1 = alpha_hat;


% Declip using ell_1 reweighted with clipping constraints
disp('Solving reweighted ell_1 ...')
w = diag(ones(N,1));
n_iters = 0;
max_iters = 10;
last_alpha_hat = zeros(N,1);
while 1,
    cvx_begin
        variable alpha_hat(N) complex;
        minimize(norm(w*alpha_hat,1));
        subject to
            Phi * Psi * alpha_hat == y;
            real(Psi_u * alpha_hat) >= cl;
            real(Psi_l * alpha_hat) <= -cl;
    cvx_end
    n_iters = n_iters + 1;
    if norm(alpha_hat - last_alpha_hat) < 0.1 || n_iters == max_iters,
        break
    end
    last_alpha_hat = alpha_hat;
    w=1./abs(alpha_hat);
    w(w>1000)=1000;
    w = diag(w);
end

x_hat_rew = real(ifft(alpha_hat)); 

err_bp = norm(x_hat_bp - x);
err_bpcc = norm(x_hat_bpcc - x);
err_rew = norm(x_hat_rew - x);

disp('Reconstruction errors')
fprintf('BP: %f, BPCC: %f, Reweighted: %f \n', err_bp, err_bpcc, err_rew)
%% Save the results
save_results = false;
if save_results,
    file_name = [mfilename '_' datestr(now,30)];
    save(file_name, '-V7')
end

%% Plot the results
clf
plot(n, x, 'k--', 'linewidth', 2)
hold on
plot(n, x_hat_bp, 'r')
plot(n, x_hat_bpcc, 'g')
plot(n, x_hat_rew, 'b')
legend('x', 'BP', 'BPCC', 'Rew')
