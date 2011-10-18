close all
clear all
clc

N = 128;
M = 60;
cvx_quiet(true)
x = get_sparse_dct_signal(N, 5);
cl = get_clip_level(x, M);
x_c = clip_signal(x, cl);   

alpha = dct(x);
alpha_c = dct(x_c);
%figure, stem(abs(alpha)), hold on, stem(abs(alpha_c), 'r')

x_hat = declip_dct(x, cl, 'consOMP');
%x_hat = declip_dct(x, cl, 'TPCC');

err = norm(x - x_hat)

figure
plot(x)
hold on
plot(x_hat, 'r')
line([1 N], [cl cl], 'color', 'g')
line([1 N], [-cl -cl], 'color', 'g')
