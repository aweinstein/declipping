% Declip an acceleration signal (from Dr. Mooney) usign TCPP

clear all
close all
clc

%%
load accel_data
offset = 1;
N = 512;
x_c = y_c(offset:offset+N);
Cu = clipu;
Cl = clipl;

N = length(x_c);
i_uc = x_c >= Cu;
i_lc = x_c <= Cl;
i_clip =  i_uc | i_lc;

% Take as measurements the values of x_c that are not clipped
I = eye(N);
Phi = I(~i_clip,:);
y = Phi * x_c;

Psi = ifft(eye(N));
A = Phi * Psi;

r = y;
alpha_c = abs(fft(x_c));
i = 1;
Delta = []; 
epsilon = 60;
%%
%figure(1)
while norm(r) > epsilon,
    k = find(alpha_c == max(alpha_c(1:N/2)), 1);
%     if new_index > 1,
%         Delta = union(Delta, [new_index, N - new_index + 2]);
%     else
%         Delta = union(Delta, new_index);
    Delta = union(Delta, [k, mod(N - k + 1, N) + 1]);
    %Delta, pause
    alpha_c(Delta) = 0;
    alpha = pinv(A(:,Delta)) * y;
    r = y - A(:,Delta) * alpha;
    i = i + 1;
%     norm(r), Delta
%     alpha_hat = zeros(N,1);
%     alpha_hat(Delta) = alpha;
%     x_hat = real(ifft(alpha_hat));
%     clf, plot(x_c), hold on, plot(x_hat,'r'), pause
end
alpha_hat = zeros(N,1);
alpha_hat(Delta) = alpha;
x_hat = real(ifft(alpha_hat));

%% Plotting
figure
plot(x_c,'b'), hold on
plot(x_hat,'r--')

save_results = false;
if save_results,
    file_name = [mfilename '_' datestr(now,30)];
    save(file_name, '-V7')
end
