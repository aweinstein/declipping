% Support estimation using the DFT
% Figure 2 of the paper

close all
clear all
clc

%%
N = 2^7;        
n = (0:N-1)';
T = N;
% M = round(0.85 * N);
M = 40;
K = 5;

paper_version = true; % Use the same signal as in the paper
if paper_version,
    amps = [1 0.65 0.3 0.825 0.475];
    fs = [17 21 42 54 62];
    x = get_sparse_signal(N,1,K,fs,amps,[]);
else
    opt.f = true; opt.phase = false; opt.A = false;
    x = get_sparse_signal(N,1,K,[],[linspace(1,0.3,K)],[], opt);
end
cl = get_clip_level(x, M);
i_uc = x >= cl;
i_lc = x <= -cl;
x_c = x;
x_c(i_uc) = cl;
x_c(i_lc) = -cl;

alpha = fft(x);
alpha_c = fft(x_c);

figure(1), clf
subplot(211)
plot(x), hold on, plot(x_c,'r')
subplot(212)
stem(abs(alpha)), hold on, stem(abs(alpha_c),'r')

save_results = false;
if save_results,
    file_name = ['mat_files/' mfilename '_' datestr(now,30)];
    save(file_name, '-V7')
end
