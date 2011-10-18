function x_c = clip_signal(x, cl)
% CLIP_SIGNAL - Clip a signal
% 
% Syntax:  x_c = clip_signal(x, cl)
%          
% Inputs:
%    x - Input signal
%    cl - clipping level
%
% Outputs:
%    x - clip signal
%    
% Example:
%  x = get_sparse_signal(N,1,2,[],[],[]);
%  x_c = clip_signal(x, 1);
%
% Other m-files required: None

% Author: Alejandro Weinstein
% Colorado School of Mines
% email: alejandro.weinstein@gmail.com
% August 2010; Last revision: 2010-08-05

cl = abs(cl);
x_c = x;
x_c(x >= cl) = cl;
x_c(x <= -cl) = -cl;
