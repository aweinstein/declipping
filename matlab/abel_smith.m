% Technique used by
%
% J.S. Abel and J.O. Smith, "Restoring a clipped signal," Acoustics, Speech, 
% and Signal Processing, 1991. ICASSP-91., 1991 International Conference on, 
% 1991, p. 1745–1748.

close all
clear all
clc

N = 2^7;        
n = (0:N-1)';
T = N;  
M = 70; % 33 is the threshold for recovery with the signal in the paper

signal_type = 1;
if signal_type == 1,
    options.f = false;
    options.A = true;
    options.phase = false;
    K = 5;
    x = get_sparse_signal(N,1,K,[],[],[],options);
    M = 40;
elseif signal_type == 2,
    load y_abel
    x = y_hat;
end
cl = get_clip_level(x, M);
clip_level = cl;
i_uc = x >= clip_level;
i_lc = x <= -clip_level;
x_c = x;
x_c(i_uc) = cl;
x_c(i_lc) = -cl;

I = eye(N);
Phi = I(~(i_uc | i_lc),:);
Phi_u = I(i_uc,:);
Phi_l = I(i_lc,:);

s = x;
c = x_c;

W = fft(eye(N));
zero_band = N/8+2:N-N/8;

% Using equiations (6) and (7)
Nc = sum(i_uc) + sum(i_lc);
Tc = zeros(N, Nc);
j = 1;
for i = 1:N, % There must be a better way to do this!
    if i_uc(i) == 1,
        Tc(:,j) = I(:,i);
        j = j + 1;
    elseif i_lc(i) == 1,
        Tc(:,j) = -I(:,i);
        j = j + 1;
    end
end

Bp = W(zero_band,:);

cvx_begin
    variable z(Nc);
    minimize(0);
    subject to
        z >= 0;
        Bp * Tc * z == -Bp * c;
cvx_end

r = c + Tc * z;

error_abel = norm(r-s)

rbptc = rank(Bp*Tc);
cond_n = cond(Bp*Tc)
fprintf('\nrank(Bp*Tc): %d Nc: %d\n', rbptc, Nc)
if rbptc >= Nc,
    disp('Theory predicts a unique solution')
else
    disp('Theory predicts multiple solutions')
end

%% Declip using TCPP
x_tpcc = declip(x, cl, 'tpcc');
err_tpcc = norm(x - x_tpcc)
%% Plot the results
plot(n,s,n,r)

line([0 N-1], [cl cl],'Color','r')
line([0 N-1], [-cl -cl],'Color','r')
ylim([min(x)*1.1 max(x)*1.1])
a_r = abs(fft(r));
figure, stem(a_r), hold on,stem(abs(fft(x)),'g')
yl = get(gca, 'YLim');
line([zero_band(1) zero_band(1)],[0 yl(2)],'color','r')
line([zero_band(end) zero_band(end)],[0 yl(2)],'color','r')
xlim([0 N+1])
legend('reconstructed', 'original')



