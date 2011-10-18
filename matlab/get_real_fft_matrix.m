function Psi = get_real_fft_matrix(N)
% GET_REAL_FFT_MATRIX - Get the corresponding FFT matrix with real entries
% 
% Syntax:  Psi = get_real_fft_matrix(N)
%
% Inputs:
%    N - Dimension of the matrix
%             
% Outputs:
%    Psi - Real FFT matrix
%
% Example: 
%    Psi = get_real_fft_matrix(128);
%
% Other m-files required: None

% Author: Alejandro Weinstein
% Colorado School of Mines
% email: alejandro.weinstein@gmail.com
% July 2010; Last revision: 2010-08-11

W = ifft(eye(N));
Psi = zeros(N,N);
Psi(:,1) = W(:,1);
Psi(:,2:N/2) = real(W(:,2:N/2));
Psi(:,N/2+1:N-1) = imag(W(:,2:N/2));
Psi(:,N) = W(:,N/2+1);
Psi = Psi * diag(1./sqrt(sum(Psi.*conj(Psi))));
