function [ q2s_mean, q2s_var ] = compute_quant_to_source( s2q_mean, s2q_var, Q, Q0, u_prob, b )
% computes the Q to S messages (cf Section 5.2.1)
%
% arguments:
%  s2q_mean:    m*n sparse matrix; mean of the S to Q Gaussian messages
%  s2q_var:     m*n sparse matrix; variance of the S to Q Gaussian messages
%  Q:           m*n sparse matrix; the quantization matrix
%  Q0:          m*1 vector; the offset vector
%  u_prob:      m*(2^b) matrix; the U to Q messages (== T to U messages), 
%               j-th row represents p(u_a == [j-th index])
%  b:           scalar; the size of the translator output
%
% returns:
%  q2s_mean:    m*n sparse matrix; mean of the Q to S Gaussian messages
%  q2s_var:     m*n sparse matrix; variance of the Q to S Gaussian messages

num_bins = 2^b;
bin = -num_bins/2:num_bins/2-1;

%% define helping variables
qpos = spones(Q); %0-1 representation of sparsity of Q
QQ = Q.^2;
Qm = Q.*s2q_mean;
Qms = sum(Qm, 2);
QQv = QQ.*s2q_var;
Qinv = spfun(@(x) 1./x, Q);
QQinv = Qinv.^2;
binmean = sum(bsxfun(@times, u_prob, bin),2);
binvar = sum(bsxfun(@times, u_prob, bin.^2),2) - binmean.^2;

%% calculate Q to S messages (Eq. 5.28 and 5.30)
q2s_mean = - Qinv .* (diag(sparse(Qms)) * qpos - Qm) + diag(sparse(- Q0 + binmean + .5)) * Qinv; 
q2s_var = QQinv .* (diag(sparse(sum(QQv, 2))) * qpos - QQv) + diag(sparse(binvar)) * QQinv + convert_slab(Qinv);

%% take care of nan
if sum(q2s_mean > 100) > 0
    error('q2s_mean > 100');
end
if sum(q2s_var > 100) > 0
    error('q2s_var > 100');
end

%% take care of inf
if sum(isinf(q2s_mean)) > 0
    q2s_mean
    error('t2s_mean inf');
end
if sum(isinf(q2s_var)) > 0
    q2s_var
    error('t2s_var inf');
end

%% take care of nan
if sum(isnan(q2s_mean)) > 0
    q2s_mean
    error('t2s_mean nan');
end
if sum(isnan(q2s_var)) > 0
    q2s_var
    error('t2s_var nan');
end

