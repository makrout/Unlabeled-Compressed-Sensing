function [MEAN,VAR, mu_u2_ext] = row_prior_U(A,B, B_u2)
    RANK=size(A,1);n=size(B,1);
    AA=repmat(diag(A),1,n);
   
    Prob=-0.5*AA+B';
    KeepMax=max(Prob);
    Prob=exp(Prob-repmat(KeepMax,RANK,1));
    
    B_u2=min(B_u2,1-1e-10);
    Prob=Prob.*B_u2'./(1-B_u2');
    
    Norm=repmat(sum(Prob),RANK,1);
    Norm = max(Norm, 1e-10);
    Prob=Prob./Norm;
    MEAN=Prob';
    VAR=diag(sum(MEAN)/n)-Prob*Prob'/n;
    
        
    % compute the extrinsic mean variable for the U+ prior
    mu_u2_ext = MEAN.*(1-B_u2)./(B_u2.*(1-MEAN)+MEAN.*(1-B_u2));

end