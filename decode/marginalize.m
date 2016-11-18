function [ s_mean, s_var ] = marginalize( s2s_prod_mean, s2s_prod_var, q2s_mean, q2s_var)
% computes the marginal of the source values (cf Section 6.2.3)
%
% arguments:
%  s2s_prod_mean:   n*1 vector; mean of product of intra-source messages
%  s2s_prod_var:	n*1 vector, variance of product of intra-source messages
%  q2s_mean:        m*n sparse matrix; mean of the Q to S Gaussian messages
%  q2s_var:         m*n sparse matrix; variance of the Q to S Gaussian messages
%
% returns:
%  s_mean:          n*1 vector: mean of the marginal (MMSE estimate)
%  s_var:           n*1 vector: variance of the marginal

varinv = spfun(@(x) 1./x, q2s_var);
s_var = 1./ (1./s2s_prod_var + sum(varinv,1)'); % Eq. 6.33

mdivv = q2s_mean.*varinv;
s_mean = (s2s_prod_mean./s2s_prod_var + sum(mdivv, 1)') .* s_var; % Eq. 6.34