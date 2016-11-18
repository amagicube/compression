function [ u_prob ] = compute_quant_to_alphabet( s2q_mean, s2q_var, Q, Q0, b )
% computes the Q to U messages (cf Section 5.2.3)
%
% arguments:
%  s2q_mean:    m*n sparse matrix; mean of the S to Q Gaussian messages
%  s2q_var:     m*n sparse matrix; variance of the S to Q Gaussian messages
%  Q:           m*n sparse matrix; the quantization matrix
%  Q0:          m*1 vector; the offset vector
%  b:           scalar; the size of the translator output
%
% returns:
%  u_prob:      m*(2^b) matrix; the Q to U messages (== U to T messages), 
%               j-th row represents p(u_a == [j-th index])

%% define helping variables
num_bins = 2^b;
m = size(Q,1);
bin_internal = -num_bins/2+1:num_bins/2-1;

QQ = Q.^2;
Qm = Q.*s2q_mean;
Qms = sum(Qm, 2);
QQv = QQ.*s2q_var;

%% calculate u_prob (Eq. 5.49)
QQvsrt = sum(QQv,2).^0.5; 
basec = (-Qms - Q0) ./ QQvsrt; % first 2 terms in Eq. 5.48
addc = bsxfun(@rdivide, bin_internal, QQvsrt); % last term in Eq. 5.48
cc = bsxfun(@plus, basec, addc); % value in Eq. 5.48
u_prob = diff([zeros(m,1), normcdf(cc), ones(m,1)], 1, 2); %Eq. 5.49

