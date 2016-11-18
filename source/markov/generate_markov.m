function [ s ] = generate_markov( model, n )
% generate the source sequence according to the Gauss-Markov model 
%  (cf Section 6.1.2)
%
% arguments:
%  model:       struct; contains information about the source
%      a:       scalar; the autocorrelation parameter (Eq. 6.5)
%      g:       scalar; variance of innovation (Eq. 6.5), = 1/(lambda^s)
%      g0:      scalar; variance of 1st element (Eq. 6.4), = 1/(lambda_0^s)
%      m0:      scalar; mean of 1st element
%  n:           scalar; the length of the sequence to be generated
%
% returns:
%  s:           n*1 vector: the generated source sequence

%% scalar version
s = zeros(n,1);
s(1) = normrnd( model.m0, sqrt(model.g0) );

for i = 2:n
    s(i) = model.a*s(i-1) + normrnd( 0, sqrt(model.g) );
end