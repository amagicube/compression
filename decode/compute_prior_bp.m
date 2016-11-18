function [ s2s_prod_mean, s2s_prod_var, s2s ] = compute_prior_bp( q2s_prod_mean, q2s_prod_var, s2s, model )
% selects the appropriate source bp decoding given the model
%  (cf Section 6.2.2)
%
% arguments:
%  q2s_prod_mean:   n*1 vector; the mean of the product of the Q to S 
%                   Gaussian messages
%  q2s_prod_var:	n*1 vector; the variance of the Q to S Gaussian
%                   messages 
%  s2s:             struct; model dependent, representing the messages
%                   within the source graph from the previous iteration
%  model:           struct; containing information about the source model
%
% returns:
%  s2s_prod_mean:   n*1 vector; the mean of the S output messages
%  s2s_prod_var:    n*1 vector; the variance of the S output messages
%  s2s:             struct; model dependent, representing the messages 
%                   within the source graph, saved for next iteration

if ~isfield(model, 'name') 
    error('model name missing');
elseif strcmp(model.name, 'iid')
    [ s2s_prod_mean, s2s_prod_var, s2s ] = compute_iid( q2s_prod_mean, q2s_prod_var, s2s, model );
elseif strcmp(model.name, 'learn')
    [ s2s_prod_mean, s2s_prod_var, s2s ] = compute_learn( q2s_prod_mean, q2s_prod_var, s2s, model );
elseif strcmp(model.name, 'markov')
    [ s2s_prod_mean, s2s_prod_var, s2s ] = forward_backward( q2s_prod_mean, q2s_prod_var, s2s, model );
else
    error('model name not recognized');
end

%% take care of nan
if sum(isnan(s2s_prod_mean)) > 0
    s2s_prod_mean
    error('s2s_prod_mean nan');
end
if sum(isnan(s2s_prod_var)) > 0
    s2s_prod_var
    error('s2s_prod_var nan');
end