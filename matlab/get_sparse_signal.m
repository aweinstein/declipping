function x = get_sparse_signal(T, n_T, K, f, A, phases, varargin)
% GET_SPARSE_SIGNAL - Generate a sparse signal in the Fourier domain
% 
% Syntax:  x = get_sparse_signal(T, n_T, K, f, A, phases)
%          x = get_sparse_signal(T, n_T, K, f, A, phases, options)
% Inputs:
%    T - Fundamental period of the signal. 
%    n_T - Number of periods to be generated.
%    K - Number of sinusoids.
%    f - Vector with the frequencies. The values are rounded to integer values 
%        between 1 and T/2, to get ``on the grid'' data.
%    A - Vector with the amplitudes of each sinusoid.
%    phases - Vector with the phases of each sinusoid, in rad.
%    options - Structure with option to randomize the parameters. The structure
%              must contain the fields options.f, options.A and options.phase.
%              Each of the fiels is boolean, where a true value indicates that 
%              the correspondig parameter is generated randomly.
% Outputs:
%    x - generated signal
%    
% Example:
%  x = get_sparse_signal(N,1,2,[],[],[]);
%
% Other m-files required: None

% Author: Alejandro Weinstein
% Colorado School of Mines
% email: alejandro.weinstein@gmail.com
% July 2010; Last revision: 2010-08-05

n_varargin = size(varargin, 2);
if  n_varargin == 0,
    opt.f = false;
    opt.A = false;
    opt.phase = false;
elseif n_varargin == 1,
    opt = varargin{1};
else
    error('Too many arguments')
end

if isempty(f),
    if opt.f,
        tmp = 2:T/2;
        rp = randperm(T/2-1);
        f(1:K) =  tmp(rp(1:K));
    else
         f = 1:K;
    end
end
if isempty(A),
    if opt.A,
        A = sign(randn(K,1)) .* (0.5+rand(K,1)); % Add 0.5 to avoid close to 0
                                                 % coefficients
        %A = 2*randn(K, 1);                                        
    else
        A = ones(K,1);
    end
end

if isempty(phases),
    if opt.phase,
        phases = 2*pi*rand(K,1);
    else
        phases = zeros(K,1);
    end
end

if length(f) ~= K || length(A) ~= K || length(phases)~=K,
    error('f, A and phases must have length K')
end

f = round(f); % round to integer values to get data in the grid
N = T * n_T;
n = (0:N-1)';
x = zeros(N,1);
for i = 1:K,
    x = x + A(i) * sin(2*pi*n*f(i)/T + phases(i));
end
