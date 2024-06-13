function [mu_u1_ext] = column_prior_U(B)
    RANK=1;
    B=max(min(B,0.1-1e-7),1e-7);
    Prob = (B./(1-B));
    Norm=repmat(sum(Prob,2),RANK,1);    
    Prob=Prob./Norm;
    MEAN=Prob';
    
    % compute the extrinsic mean variable for the U- prior
    mu_u1_ext=MEAN.*(1-B')./(B'.*(1-MEAN)+MEAN.*(1-B'));   
end