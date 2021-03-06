function x_hat = declip(x, clip_level, method)
% DECLIP - Declip a signal sparse in the frequency domain
% 
% Syntax:  x_hat = declip(x, clip_level, method)
%
% Inputs:
%    x - Original signal
%    clip_level - [cl_u cl_l] Upper and lower clipping levels. A single
%                 number indicate symmetric clipping
%    method - Method used to declip the signal. It must be one of the 
%             following:
%             'bp' : Standard Basis Pursuit
%             'bpcc' : Basis Pursuit with clipping constraints
%             'rw_ell1': Reweighted $ell_1$ minimization
%             'tcpp': Trivial Pursuit with Clipping Constraints
%
% Outputs:
%    x_hat - estimated declipped signal 
%
% Example: 
%    
%
% Other m-files required: CVX

% Author: Alejandro Weinstein
% Colorado School of Mines
% email: alejandro.weinstein@gmail.com
% July 2010; Last revision: 2011-07-13

if numel(clip_level) == 1,
    cl_u = abs(clip_level);
    cl_l = -abs(clip_level);
elseif numel(clip_level) == 2,
    cl_u = max(clip_level);
    cl_l = min(clip_level);
else
    error('clip_level must have one or two elements')
end

cvx_solver SDPT3;

switch lower(method),
    case('bp'),
        x_hat = declip_bp(x, cl_u, cl_l);
    case('bpcc'),
        x_hat = declip_bpcc(x, cl_u, cl_l);
    case('rw_ell1'),
        x_hat = declip_rw_ell1(x, cl_u, cl_l);
    case('tpcc'),
        x_hat = declip_tpcc(x, cl_u, cl_l);
    otherwise,
        error('unknown method')
end

function x_hat = declip_bp(x, cl_u, cl_l)
    N = length(x);
    [Phi, ~, ~, x_c] = get_Phi(x, cl_u, cl_l);
    y = Phi * x_c;
    Psi = get_Psi(N);
    A = Phi * Psi;
    W = diag(sqrt(sum(A.*conj(A))));
    cvx_begin
        variable alpha_hat(N) complex;
        minimize(norm(W * alpha_hat,1));
        subject to
                A * alpha_hat == y
    cvx_end
    x_hat = real(ifft(alpha_hat));

function x_hat = declip_bpcc(x, cl_u, cl_l)
    N = length(x);
    [Phi, i_uc, i_lc, x_c] = get_Phi(x, cl_u, cl_l);
    y = Phi * x_c;
    Psi = get_Psi(N);
    Psi_u = Psi(i_uc,:);
    Psi_l = Psi(i_lc,:);
    A = Phi * Psi;
    W = diag(sqrt(sum(A.*conj(A))));
    cvx_begin
        variable alpha_hat(N) complex;
        minimize(norm(W*alpha_hat,1));
        subject to
            A * alpha_hat == y;
                real(Psi_u * alpha_hat) >= cl_u;
                real(Psi_l * alpha_hat) <= cl_l;
    cvx_end
    x_hat = real(ifft(alpha_hat));
    
function x_hat = declip_rw_ell1(x, cl_u, cl_l)
    N = length(x);
    [Phi, i_uc, i_lc, x_c] = get_Phi(x, cl_u, cl_l);
    Psi = get_Psi(N);
    Psi_u = Psi(i_uc,:);
    Psi_l = Psi(i_lc,:);
    y = Phi * x_c;
    w = diag(ones(N,1));
    n_iters = 0;
    max_iters = 4;
    last_alpha_hat = zeros(N,1);
    epsilon = 1e-3;
    while 1,
        cvx_begin
            variable alpha_hat(N) complex;
            minimize(norm(w*alpha_hat,1));
            subject to
                Phi * Psi * alpha_hat == y;
                real(Psi_u * alpha_hat) >= cl_u;
                real(Psi_l * alpha_hat) <= cl_l;
        cvx_end
        n_iters = n_iters + 1;
        if norm(alpha_hat - last_alpha_hat) < 0.1 || n_iters == max_iters,
            break
        end
        if nnz(isnan(w)) || nnz(isnan(alpha_hat)),
            break;
        end
        last_alpha_hat = alpha_hat;
        w = 1 ./ (abs(alpha_hat) + epsilon);
        %w(w>1000) = 1000;
        w = diag(w);
    end

    x_hat = real(ifft(alpha_hat));

    
function x_hat = declip_tpcc(x_c, cl_u, cl_l)
    
    N = length(x_c);
    [Phi, ~, ~, x_c] = get_Phi(x_c, cl_u, cl_l);
    Psi = get_Psi(N);
    y = Phi * x_c;
    A = Phi * Psi;
    r = y;
    alpha_c = abs(fft(x_c));
    i = 1;
    Delta = [];
    epsilon = 10e-6;
    while norm(r) > epsilon,
        k = find(alpha_c == max(alpha_c(1:N/2)), 1); 
        Delta = union(Delta, [k, mod(N - k + 1, N) + 1]);
        alpha_c(Delta) = 0;
        alphac = pinv(A(:,Delta)) * y;
        r = y - A(:,Delta) * alphac;
        i = i + 1;
    end
    
    % For some low clipping values it is possible that all the non-clipped 
    % samples are 0 => r = 0 => Delta is empty
    if ~isempty(Delta), 
        alpha_hat = zeros(N,1);
        alpha_hat(Delta) = alphac;
        x_hat = real(ifft(alpha_hat));
    else
        x_hat = x_c;
    end
    
function [Phi, i_uc, i_lc, x_c] = get_Phi(x, cl_u, cl_l)
   i_uc = x >= cl_u;
   i_lc = x <= cl_l;
   i_c = i_lc | i_uc;
   I = eye(length(x));
   Phi = I(~i_c,:);
   x_c = x;
   x_c(i_uc) = cl_u;
   x_c(i_lc) = cl_l;
   
function Psi = get_Psi(N)
   Psi = ifft(eye(N));
