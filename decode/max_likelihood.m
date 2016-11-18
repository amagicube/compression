function [ u_hat ] = max_likelihood( u_prob, b )
% outputs the maximum likelihood u 
%
% arguments:
%  u_prob:      m*(2^b) matrix; the Q to U messages (== U to T messages), 
%               j-th row represents p(u_a == [j-th index])
%  b:           scalar; the size of the translator output
%
% returns:
%  u:           m*1 vector: the maximum likelihood quantized sequence u_hat

num_bins = 2^b;

[ ~, idx ] = max(u_prob,[],2);
u_hat = idx - num_bins/2 - 1;