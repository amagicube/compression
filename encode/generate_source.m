function [ s ] = generate_source( model, n )
% select the correct model to generate the source sequence
%
% arguments:
%  model:   struct; contains information about the source
%  n:       scalar; the length of the sequence to be generated
%
% returns:
%  s:       n*1 vector: the generated source sequence

if ~isfield(model, 'name') 
    error('model name missing');
elseif strcmp(model.name, 'iid')
    s = generate_iid( model, n );
elseif strcmp(model.name, 'learn')
    s = generate_learn( model, n );
elseif strcmp(model.name, 'markov')
    s = generate_markov( model, n );
else
    error('model name not recognized');
end

