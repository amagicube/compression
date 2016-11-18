function [ model ] = model_markov( a, g, g0, m0, n )
% creates the Gauss-Markov model  (cf Section 6.1.2)
%
% arguments:
%  a:           scalar; the autocorrelation parameter (Eq. 6.5)
%  g:           scalar; variance of innovation (Eq. 6.5), = 1/(lambda^s)
%  g0:          scalar; variance of 1st element (Eq. 6.4), = 1/(lambda_0^s)
%  m0:          scalar; mean of 1st element
%  n:           scalar; the length of the sequence to be generated
%
% returns:
%  model:       struct; contains information about the source
%      name:    string; name of the model
%      a:       scalar; the autocorrelation parameter (Eq. 6.5)
%      g:       scalar; variance of innovation (Eq. 6.5), = 1/(lambda^s)
%      g0:      scalar; variance of 1st element (Eq. 6.4), = 1/(lambda_0^s)
%      m0:      scalar; mean of 1st element
%      s2s:     struct; model dependent, representing the messages
%               within the source graph

model = [];
model.name = 'markov';

%% parameters
model.a = a;
model.g = g;
model.g0 = g0;
model.m0 = m0;

%% initial messages
model.s2s = [];
model.s2s.hf = zeros(n-1,1);
model.s2s.Lf = ones(n-1,1);
model.s2s.hb = zeros(n-1,1);
model.s2s.Lb = ones(n-1,1);
model.s2s.h0 = g0 * m0;
model.s2s.hn = 0;