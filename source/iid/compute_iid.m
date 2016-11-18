function [ s2s_prod_mean, s2s_prod_var, s2s ] = compute_iid( q2s_prod_mean, q2s_prod_var, s2s, model )
% computes BP of the iid model (cf Section 6.2.2)
%
% arguments:
%  q2s_prod_mean:   n*1 vector; the mean of the product of the Q to S 
%                   Gaussian messages, not used in the iid model
%  q2s_prod_var:	n*1 vector; the variance of the Q to S Gaussian
%                   messages, not used in the iid model
%  s2s:             struct; model dependent, representing the messages
%                   within the source graph from the previous iteration,
%                   not used in the iid model
%  model:           struct; containing information about the source model
%
% returns:
%  s2s_prod_mean:   n*1 vector; the mean of the S output messages
%  s2s_prod_var:    n*1 vector; the variance of the S output messages
%  s2s:             struct; not used in the iid model

n = size(q2s_prod_mean,1);

%% trivial source/prior graph
s2s_prod_mean = model.mean * ones(n,1);
s2s_prod_var = model.var * ones(n,1);
