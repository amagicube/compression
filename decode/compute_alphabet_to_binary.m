function [ z_prob ] = compute_alphabet_to_binary( u_prob, z_prob, b, binmx )
% computes the T to Z messages (cf. Section 4.4.2, Equation 4.36),
% 
% arguments:
%  u_prob:  m*(2^b) matrix; the U to T messages, 
%           j-th row represents p(u_a == [j-th index])
%  z_prob:  mb*1 vector; the Z to T messages, represents p(z_i == 1)
%  b:       scalar; the size of the translator output
%  binmx:   1 * (2^b) * b 3D-matrix; represents choice of translator
%           (cf Section 4.2)
%
% returns:
%  z_prob:  mb*1 vector: the T to Z messages, represents p(z_i == 1)

%% define helping variables
m = size(z_prob,1) / b;
ep = .0001; %small value

%% pre-processing
z_prob_reshaped = permute(reshape(z_prob, b, m), [2,3,1]) + ep; % m * 1 * b 3D-matrix [add .0001 to avoid div by 0]
u_prob_pre = abs(bsxfun(@plus, z_prob_reshaped, binmx) - 1); % m * (2^b) * b 3D-matrix

%% take product along 3rd axis (in log domain)
u_prob_pre_log = log(u_prob_pre);
logprodmx = sum(u_prob_pre_log, 3) + log(u_prob); % m * (2^b) 2D-matrix
prodmx = exp(bsxfun(@minus, logprodmx, u_prob_pre_log)); % m * (2^b) * b 3D-matrix

%% sum and normalized (Eq. 4.36)
z_prob_out_pre = sum(bsxfun(@times, prodmx, binmx), 2) ./ sum(prodmx, 2); % m * 1 * b 3D-matrix
z_prob = reshape(permute(z_prob_out_pre, [3,1,2]), [m*b, 1]); 