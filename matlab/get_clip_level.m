function cl = get_clip_level(x, M)
% GET_CLIP_LEVEL - Get the necesary clip level to get M measurements from x
% 
% Syntax:  cl = get_clip_level(x, M)
%
% Inputs:
%    x - Signal
%    M - Desired number of unclipped samples
% Outputs:
%    cl - clip level
%
% Other m-files required: None

% Author: Alejandro Weinstein
% Colorado School of Mines
% email: alejandro.weinstein@gmail.com
% July 2010; Last revision: 2010-07-26

    % Do a binary search to find the right clipping level
    max_v = max(abs(x)); 
    cl = max_v / 2; % Initial guess
    up_cl = max_v;
    low_cl = 0;
    last_Mc = NaN;
    while 1,
        Mc = sum(abs(x) < cl);
        if  Mc == M,
            break
        elseif Mc > M,
            up_cl = cl;
            cl = (cl + low_cl) / 2;
        else
            low_cl = cl;
            cl = (cl + up_cl) / 2;
        end
        if last_Mc == Mc,
            break
        end
        last_Mc = Mc;
    end
