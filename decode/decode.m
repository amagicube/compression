function [ s_hat, o ] = decode( x, dope, H, Q, Q0, model, b, bin_code, iters, convg, o )
% decodes the compressed bit stream into the reconstructed MMSE estimate
%
% arguments:
%  x:           kb*1 vector; the compressed bit stream
%  dope:        m*1 vector; represents doped bits,
%               doped bits are set to be '1' or '0' in their respective
%               locatiosn, undoped bits are set to '0.5'
%  H:           kb*mb sparse matrix; the LDPC code matrix
%  Q:           m*n sparse matrix; the quantization matrix
%  Q0:          m*1 vector; the offset vector
%  model:       struct; containing information about the source model
%  b:           scalar; the size of the translator output
%  bin_code:    string; the translator choice: 'bin' for standard binary
%               code, 'pam' for PAM gray code, 'qam' for QAM gray code
%  iters:       scalar; maximum number of iterations to run BP
%  convg:       scalar: convergence criterion, terminate if successive
%               s_hat satisfies Eq. 6.35, ie.
%                   || s^(\tau) - s^(\tau-1) ||_\infty < convg
%  o:           struct; to pass additional information into the function
% 
% returns:
%  s_hat:       n*1 vector; the decoded value
%  o:           struct: to return additional information from the function

%% define helper variables
[m,n] = size(Q);
kb = size(x,1);

num_bins = 2^b;
if strcmp(bin_code, 'bin')
    binary_map = 0:num_bins-1;
else
    binary_map = bin2gray(0:num_bins-1,bin_code,num_bins);
end
binmx = reshape(de2bi(binary_map,b,'left-msb'), [1,num_bins,b]); % 1 * (2^b) * b 3D-matrix

%% initialize messages
if ~isfield(o,'isset')
    o.isset = 1;
end

s2s = model.s2s;
s2q_mean = spalloc(m,n,nnz(Q));
s2q_var = spones(Q);
z_prob = 0.5*ones(m*b,1);
h2z = spalloc(kb,m*b,nnz(H));

s_hat_prev = zeros(n,1);

%% Run Loopy BP for iters number of iteration, or till convergence
for i=1:iters
    [ u_prob ] = compute_quant_to_alphabet( s2q_mean, s2q_var, Q, Q0, b );
    [ z_prob ] = compute_alphabet_to_binary( u_prob, z_prob, b, binmx );
    [ z_prob, h2z ] = compute_code_bp( dope, x, H, z_prob, h2z );
    [ u_prob ] = compute_binary_to_alphabet( z_prob, b, binmx );
    [ q2s_mean, q2s_var ] = compute_quant_to_source( s2q_mean, s2q_var, Q, Q0, u_prob, b );
    [ q2s_prod_mean, q2s_prod_var ] = compute_source_to_prior( q2s_mean, q2s_var );
    [ s2s_prod_mean, s2s_prod_var, s2s ] = compute_prior_bp( q2s_prod_mean, q2s_prod_var, s2s, model );
    [ s2q_mean, s2q_var ] = compute_source_to_quant( s2s_prod_mean, s2s_prod_var, q2s_mean, q2s_var, Q );

    %% check convergence 
    [ s_hat, ~ ] = marginalize( s2s_prod_mean, s2s_prod_var, q2s_mean, q2s_var);
    s_diff = max(abs(s_hat_prev - s_hat));
    if isfield(o,'verbose') && o.verbose
        disp([num2str(i) ': ' num2str(s_diff)]);
    end
    if s_diff < convg
        break
    end
    s_hat_prev = s_hat;
end

o.s_diff = s_diff;
o.convg_iter = i;

%% take care of nan
if sum(isnan(s_hat)) > 0
    s_hat
    error('s_mean nan');
end
