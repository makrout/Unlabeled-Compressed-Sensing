function [MEAN,VAR] = prior_gauss_X(A,B,mu_gauss, var_gauss)

   gamma_gauss = mean(1./var_gauss);
    
   % clip the eigen values of A before inversing it in the svd form
   [U,D] = eig(A);
   gamma_gauss_positive = abs(diag(D))  + gamma_gauss;
   VAR = U * (kron(1./gamma_gauss_positive, ones(1,size(U,1))) .* U');
   MEAN=(B+mu_gauss*gamma_gauss)*VAR;
   VAR=VAR - diag(diag(VAR))+abs(diag(diag(VAR)));
   if min(diag(VAR))<0 
       error('The covariance matrix has at least one negative eigen value');
   end   
end
