function [ q2s_prod_mean, q2s_prod_var ] = compute_source_to_prior( q2s_mean, q2s_var )
% takes the product of the Q to S messages (cf Section 6.2.2)
%
% arguments:
%  q2s_mean:        m*n sparse matrix, mean of the Q to S Gaussian messages
%  q2s_var:         m*n sparse matrix, variance of the Q to S Gaussian messages
%
% returns:
%  q2s_prod_mean:   n*1 vector, mean of the product of the Q to S messages
%  q2s_prod_var:	n*1 vector, variance of the product of the Q to S messages

varinv = spfun(@(x) 1./x, q2s_var);
q2s_prod_var = full(1./ sum(varinv,1)');

mdivv = q2s_mean.*varinv;
q2s_prod_mean = full(sum(mdivv, 1)') .* q2s_prod_var;


%% take care of nan
if sum(isnan(q2s_prod_mean)) > 0
    q2s_prod_mean
    full(sum(mdivv, 1)')
    error('s2c_mean nan');
end
if sum(isnan(q2s_prod_var)) > 0
    q2s_prod_var
    error('s2c_var nan');
end