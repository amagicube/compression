function [ z_prob, h2z ] = compute_code_bp( dope, x, H, z_prob, h2z )
% computes message passing inside the code graph (cf Section 4.4.1)
%
% arguments:
%  dope:    m*1 vector; represents doped bits,
%           doped bits are set to be '1' or '0' in their respective
%           locatiosn, undoped bits are set to '0.5'
%  x:       kb*1 vector; the compressed bit stream
%  H:       kb*mb sparse matrix; the LDPC code matrix
%  z_prob:  mb*1 vector; the T to Z messages, represents p(z_i == 1)
%  h2z:     kb*mb sparse matrix; represents H to Z message from last iter
%
% returns:
%  z_prob:  mb*1 vector; the Z to T messages, represents p(z_i == 1)
%  h2z:     kb*mb sparse matrix; H to Z message for next iter

%% define helper variables
infcap = 100; %large value
ep = .0001; %small value

%% read input
dope_msg = log(1-dope) - log(dope); 
phi = log(1-z_prob) - log(z_prob); % Eq. 4.11
prior = (dope_msg+phi)' + sum(h2z, 1); 

%% compute z2h (Eq. 4.12)
z2h = H * diag(sparse(min(max(prior,-infcap),infcap))) - h2z + ep*H; %cap off at -10,10 to avoid inf, add .0001 to avoid div by 0 later
if sum(isnan(z2h(:)))
    error('z2x is nan');
end

%% compute h2z (Eq. 4.25)
tanhm = tanh(z2h/2);
[rowIdx,~,values] = find(tanhm);
prodtanhm = accumarray(rowIdx,values,[],@prod);

prodtanhmdivt = diag(sparse(prodtanhm)) * spfun(@(x) 1./x, tanhm);
multbyx = 2*atanh(diag(sparse(1-2*x)) * prodtanhmdivt); % Eq. 4.25

h2z = min(max(multbyx ,-infcap),infcap); %cap off at -10,10 to avoid inf,

if sum(isnan(h2z(:)))
    error('x2z is nan');
end

%% out messages (Eq. 4.33)
z_prob = (1 - tanh((sum(h2z, 1)' + dope_msg)/2)) / 2;
if sum(isnan(z_prob(:)))
    error('z_prob is nan');
end

