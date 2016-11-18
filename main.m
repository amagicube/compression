startup;

%% parameters
n = 500; % length of original data sequence
m = n;
d = '0.4x3/0.6x2'; % degree distribution

b = 4; % quantization parameter: 2^b quantization levels (Section 3.1.2)
w = .6; % width of each quantizer (Section 5.1.2.2)

total_rate = 2.7; % final compression rate (bits per sample, Equation 3.37)

rand_dope = 0; % 0 for lattice doping, 1 for random doping
rand_dope_rate = .34; % random dope rate (Section 4.3.1)
lat_dope_bits = [4]; % lattice doping (Section 4.3.3 and 4.3.4): bits to be doped, eg. [1 4] for multiple lattice doping of bits 1 and 4

bin_code = 'pam'; % choice of translator (Section 4.2), 'bin' or 'pam'

verbose = 0; % 1 for debug prints
max_iter = 150; % max number of iterations (Section 6.2.3)
convg_epsilon = .01; % convergence criterion (Section 6.2.3)

savepath = ''; % file path to save computation history as a .mat file

%% source model
model = model_markov(.7, .51, 1, 0, n); % Gauss Markov with parameters (a,g,g0,m0) of length n
% model = model_iid(0, 1);

%% generate quantizer matrix
q_option = [];
q_option.w = w;
[ Q, Q0 ] = generate_q( n, m, q_option );

%% generate LDPC
mb = m*b;
if rand_dope
    num_dope_bits = ceil(rand_dope_rate * mb);
    dope_indices = sort(randsample(mb, num_dope_bits))';   % random doping
else
    dope_indices = reshape(bsxfun(@plus, 0:b:(mb-1), lat_dope_bits'), 1,[]);
    num_dope_bits = size(dope_indices,2);
end
kb = total_rate*n - num_dope_bits;
H = generate_ldpc_wrapper( kb, mb, d, [], ldpc_option );

%% Generate source sequence
s = generate_source( model, n );
u = quantize_slice( s, Q, Q0, b );
z = translate( u, b, bin_code );
x = hash_ldpc( z, H );

%% Doping
dope = 0.5*ones(mb,1);
dope(dope_indices,1) = (z(dope_indices)==1);

%% decode
o = [];
o.verbose = verbose;

tic;
[ s_hat, o ] = decode( x, dope, H, Q, Q0, model, b, bin_code, max_iter, convg_epsilon, o );

%% print results
s_convg_iter = o.convg_iter;
s_mse = sum((s-s_hat).^2)/size(s,1);
sqnr = -10*log10(s_mse);

u_hat = quantize_slice( s_hat, Q, Q0, b );
z_hat = translate( u_hat, b, bin_code );
x_hat = hash_ldpc( z_hat, H );

z_err = nnz(z - z_hat) / mb;

disp(['b=' num2str(b) ', w=' num2str(w) ', d=' num2str(num_dope_bits/mb) ', k=' num2str(kb/mb) ', (R=' num2str((num_dope_bits+kb)/n) '), (SQNR=' num2str(sqnr) '), z_err=' num2str(z_err) ', i=' num2str(s_convg_iter)]);

toc