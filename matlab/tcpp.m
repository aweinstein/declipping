% Trivial Pursuit with Clipping Constraints 
% Figure 3 of the paper with 
%     signal_type = 1;
%     add_noise = false;
% Figure 6 of the paper with ...
close all
clear all
clc

%%
N = 2^7;
n = (0:N-1)';
T = N;

signal_type = 2;
if signal_type == 1,
    x = sin(2*pi/T*n) + 0.25*sin(2*pi/T*3*n);
    clip_level = 0.2; % We can go down to 0.2 and still recover x
elseif signal_type == 2,
    K = 10;
    opt.f = true; opt.phase = true; opt.A = true;
    x = get_sparse_signal(N,1,K,[],[],[], opt);
    M = 50;
    clip_level = get_clip_level(x, M);
elseif signal_type == 3,
    load y_abel
    x = y_hat;
    M = 68; % We can go down to 68 and still recover x
    clip_level = get_clip_level(x, M);
else
    Omega = [3 7 11 51 55 59 63 115 119 123];
    Omega = 65 + (Omega -3)/2;
    alpha = zeros(N,1);
    alpha(Omega) =  0.25 + 2*rand(numel(Omega),1) ;
    x = Psi * alpha;
end

add_noise = false;
if add_noise,
    sigma = 0.05;
    x_original = x;
    z = sigma * randn(size(x));
    x = x + z;
    epsilon = 8 * sigma;
    %epsilon = norm(z);
else
    epsilon = 1e-6;
end

i_uc = x >= clip_level;
i_lc = x <= -clip_level;
i_clip =  i_uc | i_lc;
x_c = x;
x_c(i_clip) = sign(x_c(i_clip)) .* clip_level;  

% Take as measurements the values of x that are not clipped
I = eye(N);
Phi = I(~i_clip,:);
y = Phi * x_c;
Psi = ifft(eye(N));
A = Phi * Psi;

r = y;
alpha_c = abs(fft(x_c));
i = 1;
Delta = []; 
while norm(r) > epsilon,
    k = find(alpha_c == max(alpha_c(1:N/2)), 1); 
    Delta = union(Delta, [k, mod(N - k + 1, N) + 1]);
    alpha_c(Delta) = 0;
    alpha = pinv(A(:,Delta)) * y;
    r = y - A(:,Delta) * alpha;
    i = i + 1;
end
Delta;
alpha_hat = zeros(N,1);
alpha_hat(Delta) = alpha;
x_hat = real(ifft(alpha_hat));

error_tpcc = norm(x - x_hat)
peak = max(abs(x))
peak_over_norm = max(abs(x)) / norm(x)

%% Plotting
if add_noise,
    figure, stem(abs(Psi\x)), hold on, stem(abs(Psi\x_c),'r')
    stem(abs(Psi\x_original),'g')
    legend('Noisy', 'Clipped', 'Original')
else
    figure, stem(abs(Psi\x)), hold on, stem(abs(Psi\x_c),'r') 
end
figure,plot(n,x,'b'), hold on,plot(n,x_hat,'r--'), plot(n,x_c,'g')
legend('Original','Recovered','Clipped')
title('Iterative method')

save_results = false;
if save_results,
    file_name = ['mat_files/' mfilename '_' datestr(now,30)];
    save(file_name, '-V7')
end
