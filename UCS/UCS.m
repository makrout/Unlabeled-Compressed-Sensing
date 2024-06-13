 function [u, x, nrmses] = UCS(Y, U, X, A, var_w, opt)
% UCS as a factorization Belief-Propagation solver for the unlabeled compressed sensing.
% Given the observation model Y = UAX + W, and the matrices Y and A, find:
%  - the permutation matrix U
%  - the matrix X

% Inputs :
% Y                     NxM observation matrix
% U                     the true permutation U matrix of size NxN, only used to compute the error
% X                     the true X matrix of size MxR, only used to compute the error
% A                     the sensing matrix of size MxR
% var_w                 noise variance
% opt                   object containing the experiment parameters (see the UCS_opt object)
%
% Outputs :
% u                     final estimate of the U matrix
% x                     final signal estimate of the X matrix
% nrmses                vector of NRMSEs over iterations
          
    [n,m] = size(Y);
    r = size(X,2);
    
    % initialize U and X variables
    u = randn(n,n);
    mu_Uminus_ext = 0.5 * ones(n,n);
    x=randn(m,r);
    
    % initialize the estimation variables u, v and x
    u_old = zeros(n,n);
    v_old = zeros(m,n);
    u_var = zeros(n,n);
    x_var = zeros(r,r);
    
    % initialize the bi-LMMSE variables
    A_u = zeros(n,n);
    B_u = zeros(n,n);
    A_x = zeros(r,r);
    B_x = zeros(m,r);
    mu_Uminus_ext_old = 0.5 * ones(n,n);
    
    % iteration variables
    var_Y = sum(Y.^2, 'all')/prod(size(Y));
    damp=opt.damping;
    keep_iter=true;
    t=0;
    nrmses = [];    
    
    while (t<opt.nb_iter)
        %% bi-LMMSE module for U- and X-
        % save old bi-LMMSE variables
        A_u_old = A_u;      A_x_old = A_x;
        B_u_old = B_u;      B_x_old = B_x;      
        
        % change of variable from X to V
        v = (A * x')';
        v_var = A * x_var * A';
        
        % update the bi-LMMSE variables
        B_u_new = (Y*v)/var_w - m*(u_old*v_var)*var_Y/var_w^2;
        A_u_new = (v'*v)/var_w + m*v_var/(opt.beta*var_w) - m*v_var*var_Y/var_w^2;

        B_v_new = (Y'*u)/var_w - n*(v_old*u_var)*var_Y/var_w^2;
        A_v_new = (u'*u)/var_w + n*u_var/(opt.beta*var_w) - n*u_var*var_Y/var_w^2;
        
        % change of variable from V to X
        B_x_new = B_v_new * A;
        A_x_new = A' * A_v_new * A;
        
        % save old estimates for damping and to check convergence
        u_old = u;
        v_old = v;

        % damp the bi-LMMSE variables
        A_u = (1-damp)*A_u_old + damp*A_u_new;
        A_x = (1-damp)*A_x_old + damp*A_x_new;
        B_u = (1-damp)*B_u_old + damp*B_u_new;
        B_x = (1-damp)*B_x_old + damp*B_x_new;
        
        %% the assignment prior on U- (inside the bi-LMMSE module)
        % run the row-wise assignment prior
        [u, u_var, mu_Uplus_ext] = row_prior_U(A_u, B_u, mu_Uminus_ext');
        
        %% The column-prior on U+
        % run the column-wise prior
        [mu_Uminus_ext] = column_prior_U(mu_Uplus_ext');
        
        % damp the message from the column-wise prior the row-wise prior
        mu_Uminus_ext = (1-damp)*mu_Uminus_ext_old + damp*mu_Uminus_ext;
        
        % save old estimate
        mu_Uminus_ext_old = mu_Uminus_ext;

        %% The i.i.d. prior on X+
        % compute the message from the prior module of X
        [x, x_var] = prior_gauss_X(A_x,B_x, zeros(m,r), 1);
        
        %% check convergence and compute metrics
        % check the convergence
        diff = mean(abs(v-v_old), 'all') + mean(abs(u-u_old), 'all');
        
        % compute the nrmse at iteration t
        nrmse = sqrt(mean((U*A*X' - u*A*x').^2/mean((U*A*X').^2), 'all'));
        nrmses = [nrmses nrmse];
        
        if mod(t,100)==0
            fprintf(1,'[t=%d] nrmse = %f \n', t, nrmse);
        end
        t=t+1;
        
        if (diff == 0) && (keep_iter)
            opt.nb_iter = t + 1000;
            keep_iter = false;
        end
        
    end
end


