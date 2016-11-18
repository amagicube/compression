function [ z_prob, h2z, b2x ] = compute_bsc_bp( dope, x, H, z_prob, h2z, y, B, b2x )
% computes message passing of the joint code-BSC graph (cf Section 9.4.2.2)
%
% arguments:
%  dope:    m*1 vector; represents doped bits,
%           doped bits are set to be '1' or '0' in their respective
%           locatiosn, undoped bits are set to '0.5'
%  x:       kb*1 vector; the compressed bit stream
%  H:       kb*mb sparse matrix; the LDPC code matrix
%  z_prob:  mb*1 vector; the T to Z messages, represents p(z_i == 1)
%  h2z:     kb*mb sparse matrix; represents H to Z message from last iter
%  y:       lkb*1 vector; the parity check bits on x
%  B:       lkb*kb sparse matrix; the LDGM parity check matrix
%  b2x:     lkb*kb sparse matrix; represents B to X message from last iter
%
% returns:
%  z_prob:  mb*1 vector; the Z to T messages, represents p(z_i == 1)
%  h2z:     kb*mb sparse matrix; H to Z message for next iter
%  b2x:     lkb*kb sparse matrix; B to X message for next iter

%% define helper variables
infcap = 100;
ep = .0001;

%% read input
dope_msg = log(1-dope) - log(dope); 
phi = log(1-z_prob) - log(z_prob); % Eq. 4.11
prior = (dope_msg+phi)' + sum(h2z, 1) ; 

xin = log(1-x) - log(x);

%% compute z2h (Eq. 4.12)
z2h = H * diag(sparse(min(max(prior,-infcap),infcap))) - h2z + ep*H; %cap off at -10,10 to avoid inf, add .0001 to avoid div by 0 later

%% compute x2b
tanhm = tanh(z2h/2);
[rowIdx,~,values] = find(tanhm);
prodtanhm = accumarray(rowIdx,values,[],@prod);

xdown = 2*atanh(prodtanhm);
accdown = xdown' + sum(b2x, 1) + xin'; 

x2b = bsxfun(@times, B, min(max(accdown,-infcap),infcap)) - b2x + ep*B; %cap off at -10,10 to avoid inf, add .0001 to avoid div by 0 later

tanhmxy = tanh(x2b/2);
[rowIdx,~,values] = find(tanhmxy);
prodtanhmxy = accumarray(rowIdx,values,[],@prod);
prodtanhmdivtxy = diag(sparse(prodtanhmxy)) * spfun(@(x) 1./x, tanhmxy);
twoatanhxy = 2*atanh(diag(sparse(1-2*y)) * prodtanhmdivtxy);
b2x = min(max(twoatanhxy ,-infcap),infcap); %cap off at -10,10 to avoid inf,

xup = sum(b2x, 1)' + xin;

%% compute h2z (Eq. 4.25)
prodtanhmdivt = diag(sparse(prodtanhm .* tanh(xup/2))) * spfun(@(x) 1./x, tanhm);
twoatanh = 2*atanh(prodtanhmdivt);
h2z = min(max(twoatanh ,-infcap),infcap); %cap off at -10,10 to avoid inf,


%% out messages (Eq. 4.33)
z_prob = (1 - tanh((sum(h2z, 1)' + dope_msg)/2)) / 2;


