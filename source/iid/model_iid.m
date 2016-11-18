function [ model ] = model_iid( iid_mean, iid_var )
% creates the iid model  (cf Section 6.1.1)
%
% arguments:
%  iid_mean:    scalar; mean of each element
%  iid_var:     scalar; variance of each element
%
% returns:
%  model:       struct; contains information about the source
%      name:    string; name of the model
%      mean:    scalar; mean of each element
%      var:     scalar; variance of each element
%      s2s:     struct; model dependent, representing the messages
%               within the source graph, not used in the iid model

model = [];
model.name = 'iid';

%% parameters
model.mean = iid_mean;
model.var = iid_var;

%% initialize messages
model.s2s = [];

