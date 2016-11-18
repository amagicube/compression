function [ s2q_mean, s2q_var ] = compute_source_to_quant( s2s_prod_mean, s2s_prod_var, q2s_mean, q2s_var, Q )
% computes the S to Q messages (cf Section 5.2.2)
%
% arguments:
%  s2s_prod_mean:   n*1 vector; mean of product of intra-source messages
%  s2s_prod_var:	n*1 vector; variance of product of intra-source messages
%  q2s_mean:        m*n sparse matrix; mean of Q to S Gaussian messages
%  q2s_var:         m*n sparse matrix; variance of Q to S Gaussian messages
%  Q:               m*n sparse matrix; the quantization matrix
%
% reutrns:
%  s2q_mean:        m*n sparse matrix, mean of the S to Q Gaussian messages
%  s2q_var:         m*n sparse matrix, variance of the S to Q Gaussian messages

qpos = spones(Q); %0-1 representation of sparsity of Q

varinv = spfun(@(x) 1./x, q2s_var);
s2q_var = spfun(@(x) 1./x, bsxfun(@times, qpos, 1./s2s_prod_var' + sum(varinv,1)) - varinv);

mdivv = q2s_mean.*varinv;
s2q_mean = (bsxfun(@times, qpos, (s2s_prod_mean./s2s_prod_var)' + sum(mdivv, 1)) - mdivv) .* s2q_var;


%% take care of nan
if sum(isnan(s2q_mean)) > 0
    s2q_mean
    error('s2t_mean nan');
end
if sum(isnan(s2q_var)) > 0
    s2q_var
    error('s2t_var nan');
end