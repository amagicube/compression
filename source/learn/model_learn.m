function [ model ] = model_learn( L_iid, h_theta, L_theta )
% creates the iid model with unknown theta (cf Section 9.3.1)
%
% arguments:
%  L_iid:       scalar; precision (lambda) of each element
%  h_theta:     scalar; potential (eta) of the unknown potential parameter 
%               theta
%  L_theta:     scalar; precision (lambda) of the unknown potential 
%
% returns:
%  model:       struct; contains information about the source
%      name:    string; name of the model
%      L_iid:   scalar; precision (lambda) of each element
%      theta:   scalar; potential (eta) of each element, parametrized as
%               theta
%      s2s:     struct; model dependent, representing the messages
%               within the source graph, not used in the iid-learn model

model = [];
model.name = 'learn';

%% parameters

theta_var = 1 / L_theta;
theta_mean = h_theta / L_theta;

model.theta = normrnd(theta_mean, sqrt(theta_var));
model.L_iid = L_iid;

%% initialize messages
model.s2s = [];

