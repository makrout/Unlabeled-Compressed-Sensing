% UCS_opt  : Function to set UCS_opt to their default values
%
%    Details of the options:
%   .nb_iter            max number of iterations
%   .damping            damping coefficient
%   .beta               temperature parameter to control the problem type
%                       i.e. estimation with beta = 1 and MAP with beta = inf

function opt = UCS_opt()
       opt.nb_iter=100000;
       opt.damping=0.2;
       opt.beta = 1;
end