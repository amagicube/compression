function [ s2s_prod_mean, s2s_prod_var, s2s ] = compute_learn( q2s_prod_mean, q2s_prod_var, s2s, model )
% computes BP of the iid model augmented by a parameter learning subgraph
%  (cf Section 9.3.1)
%
% arguments:
%  q2s_prod_mean:   n*1 vector; the mean of the product of the Q to S 
%                   Gaussian messages
%  q2s_prod_var:	n*1 vector; the variance of the Q to S Gaussian
%                   messages 
%  s2s:             struct; model dependent, representing the messages
%                   within the source graph from the previous iteration,
%                   not used in the iid-learn model
%  model:           struct; containing information about the source model
%
% returns:
%  s2s_prod_mean:   n*1 vector; the mean of the S output messages
%  s2s_prod_var:    n*1 vector; the variance of the S output messages
%  s2s:             struct; model dependent, representing the messages 
%                   within the source graph, saved for next iteration,
%      theta_mean:  scalar; mean of marginal of estimated parameter theta
%      theta_var:   scalar; var of marginal of estimated parameter theta

%% helper variables
n = size(q2s_prod_mean,1);

%% potentials
Lin = 1./q2s_prod_var;
hin = q2s_prod_mean .* Lin;

Li = model.L_iid;
hi = 0;
Lt = model.L_theta + n/model.L_iid;
ht = model.h_theta;
Lit = -1;

%% parallel messages
hu = -Lit .* (hi + hin) ./ (Li + Lin);
Lu = -Lit^2 ./ (Li + Lin);

hu_sum = sum(hu,1);
Lu_sum = sum(Lu,1);

hd = -Lit .* (ht + hu_sum - hi) ./ (Lt + Lu_sum - Li);
Ld = -Lit^2 ./ (Li + Lu_sum - Li);

%%
%% marginalize
Li_hat = Li + Ld;
hi_hat = hi + hd;

Lt_hat = Lt + Lu_sum;
ht_hat = ht + hu_sum;

%% compute output
s2s_prod_var = 1./Li_hat;
s2s_prod_mean = s2s_prod_var .* hi_hat;

%% save messages
s2s.theta_mean = ht_hat / Lt_hat;
s2s.theta_var = 1 / Lt_hat;

