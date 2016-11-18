function [ u_prob ] = compute_binary_to_alphabet( z_prob, b, binmx )
% computss the T to U messages (cf Section 4.4.2, Equation 4.35)
%
% arguments:
%  z_prob:  mb*1 vector; the Z to T messages, represents p(z_i == 1)
%  b:       scalar; the size of the translator output
%  binmx:   1 * (2^b) * b 3D-matrix; represents choice of translator
%           (cf Section 4.2)
% returns:
%  u_prob:  m*(2^b) matrix; the T to U messages, 
%           j-th row represents p(u_a == [j-th index])

%% define helping variables
m = size(z_prob,1) / b;

%% pre-processing
z_prob_reshaped = permute(reshape(z_prob, b, m), [2,3,1]); % m * 1 * b 3D-matrix
u_prob_pre = abs(bsxfun(@plus, z_prob_reshaped, binmx) - 1); % m * (2^b) * b 3D-matrix

%% take product along 3rd axis (Eq. 4.35)
u_prob = exp(sum(log(u_prob_pre), 3)); 