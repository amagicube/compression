function [ s2s_prod_mean, s2s_prod_var, s2s ] = forward_backward( q2s_prod_mean, q2s_prod_var, s2s, model )
% computes BP of the Gauss-Markov model, ie. the forward-backward algorithm
%  (cf Section 6.2.2)
%
% arguments:
%  q2s_prod_mean:   n*1 vector; the mean of the product of the Q to S 
%                   Gaussian messages
%  q2s_prod_var:	n*1 vector; the variance of the Q to S Gaussian
%                   messages 
%  s2s:             struct; model dependent, representing the messages
%                   within the source graph from the previous iteration,
%                   for Gauss-Markov it is the forward backward messages
%      Lf:          (n-1)*1 vector; the precision parameter (lambda) of the
%                   forward messages (Eq. 6.23)
%      Lb:          (n-1)*1 vector; the precision parameter (lambda) of the
%                   backward messages (Eq. 6.23)
%      hf:          (n-1)*1 vector; the potential parameter (eta) of the
%                   forward messages (Eq. 6.23)
%      hb:          (n-1)*1 vector; the potential parameter (eta) of the
%                   backward messages (Eq. 6.23)
%      h0:          scalar; the node potential (eta) of the first node 
%                   (Eq. 6.19)
%      hn:          scalar; the node potential (eta) of the last node 
%                   (Eq. 6.19)
%  model:           struct; containing information about the source model
%
% returns:
%  s2s_prod_mean:   n*1 vector; the mean of the S output messages
%  s2s_prod_var:    n*1 vector; the variance of the S output messages
%  s2s:             struct; model dependent, representing the messages 
%                   within the source graph, saved for next iteration,
%                   for Gauss-Markov it is  the forward backward messages
%      Lf:          (n-1)*1 vector; the precision parameter (lambda) of the
%                   forward messages (Eq. 6.23)
%      Lb:          (n-1)*1 vector; the precision parameter (lambda) of the
%                   backward messages (Eq. 6.23)
%      hf:          (n-1)*1 vector; the potential parameter (eta) of the
%                   forward messages (Eq. 6.23)
%      hb:          (n-1)*1 vector; the potential parameter (eta) of the
%                   backward messages (Eq. 6.23)

%% helper variables
n = size(q2s_prod_mean,1);
a = model.a; % nn*nn matrix
g = model.g; % nn*nn matrix
g0 = model.g0; % nn*nn matrix

Lf = s2s.Lf;
Lb = s2s.Lb;
hf = s2s.hf;
hb = s2s.hb;

h0 = s2s.h0;
hn = s2s.hn;

%% potentials
Lin = 1./q2s_prod_var;
hin = q2s_prod_mean .* Lin;

Li = [1/g0; 1/g*ones(n-1,1)] + [a^2/g*ones(n-1,1); 0];
hi = [h0; zeros(n-2,1); hn];
Lij = -a/g;

Lpass = Lin + Li;
hpass = hin + hi;

%% parallel messages
hf = -Lij ./ (Lpass(1:n-1) + [0;Lf(1:n-2)]) .* (hpass(1:n-1) + [0;hf(1:n-2)]);
Lf = -Lij.^2 ./ (Lpass(1:n-1) + [0;Lf(1:n-2)]);

hb = -Lij ./ (Lpass(2:n) + [Lb(2:n-1);0]) .* (hpass(2:n) + [hb(2:n-1);0]);
Lb = -Lij.^2 ./ (Lpass(2:n) + [Lb(2:n-1);0]);

%% marginalize
Li_hat = Li + [0;Lf] + [Lb;0];
hi_hat = hi + [0;hf] + [hb;0];

%% compute output
s2s_prod_var = 1./Li_hat;
s2s_prod_mean = s2s_prod_var .* hi_hat;

%% save messages
s2s.Lf = Lf;
s2s.Lb = Lb;
s2s.hf = hf;
s2s.hb = hb;