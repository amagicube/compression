function [ s ] = generate_learn( model, n )
% generate the source sequence with the iid model with unknown theta
%  (cf Section 9.3.1)
%
% arguments:
%  model:       struct; contains information about the source
%      L_iid:   scalar; precision (lambda) of each element
%      theta:   scalar; potential (eta) of each element, parametrized as
%               theta
%  n:           scalar; the length of the sequence to be generated
%
% returns:
%  s:           n*1 vector: the generated source sequence

s_var = 1 / model.L_iid;
s_mean = model.theta / model.L_iid;

s = normrnd( s_mean , sqrt(s_var), [n,1] );

