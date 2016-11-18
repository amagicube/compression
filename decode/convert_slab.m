function [ out_var ] = convert_slab( qwidth )
% computes the variance adjustment given the width, which is equal to the
%  variance of the uniform distribution with the same width 
%  (cf Equation 5.22)
%
% arguments:
%  qwidth:  scalar; the width of the slab
%
% returns:
%  out_var: scalar; the variance adjustment

out_var = qwidth.^2 / 12;