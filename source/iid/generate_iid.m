function [ s ] = generate_iid( model, n )
% generate the source sequence according to the iid model (Section 6.1.1)
%
% arguments:
%  model:       struct; contains information about the source
%      mean:    scalar; mean of the iid distribution
%      var:     scalar; var of the iid distribution
%  n:           scalar; the length of the sequence to be generated
%
% returns:
%  s:           n*1 vector: the generated source sequence

s = normrnd( model.mean, sqrt(model.var), [n,1] );

