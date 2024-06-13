close all;
clear all;
clear variables;
clc;

addpath('priors/');
addpath('UCS/');

% Problem size
n=50;
m=100;
r=10;
fprintf(1,'===== Problem dimension =====\n - N = %d\n - M = %d\n - R = %d\n\n',m,n,r);

% SNR level in dB
SNR_dB = 30;

% Variance annealing factor
var_annealing = 20;

% create the options object
opt = UCS_opt();

% define the permutation matrix
U = eye(n);
U = U(randperm(n),:);
% define the matrix X
X = randn(m,r);
% define the sensing matrix
A = randn(n,r);   

% Define the Wless Z
Z = U*A*X';

% Define the noise variance corresponding the SNR level
var_Z = sum(Z.^2, 'all')/prod(size(Z));
var_w = var_Z * 10^(-SNR_dB/10);

% Define the W vector
W = sqrt(var_w)*randn(n,m);

Y =  Z + W;

% Run BiVamp for UCS
disp('===== Running UCS =====')

tstart = tic;
[u_est, v_est, nrmses] = UCS(Y, U, X, A, var_w*var_annealing, opt);
tUCS = toc(tstart);

fprintf(1,'Final nrmse = %f \n', nrmses(end));
fprintf(1,'Running time = %f \n', tUCS);

%% Plotting
f = figure;
f.Position = [100 100 800 700];

% plot the U matrix and its estimate 
subplot(2,2,1)
[~, max_index] = max(u_est);
u_max = zeros(n);
for ii=1:n
    u_max(max_index(ii), ii) = 1;
end
stem(U(:),'b');
hold on;
stem(u_max(:),'r','--');
legend('true', 'estimated')
title({'True and estimated permutation','matrices $\mathbf{U}$ and $\widehat{\mathbf{U}}$'}, 'interpreter','latex');

% plot the X matrix and its estimate
subplot(2,2,2)
stem(X(:,1),'b');
hold on;
stem(v_est(:,1),'r','--');
legend('true', 'estimated');
title({'True and estimated', 'matrices $\mathbf{X}$ and $\widehat{\mathbf{X}}$'}, 'interpreter','latex');

% the NRMSE plot
subplot(2,2,[3,4]);
semilogy(1:length(nrmses), nrmses);
xlabel('iterations');
ylabel('NRMSE');
grid on;
title('NRMSE over iterations');
