function [ H ] = generate_ldpc_wrapper( kb, mb, d, omit_hash_indices, ldpc_option )
% MATLAB wrapper for calling MEX file
% 
% arguments:
%  kb:                  scalar; the height of the LDPC matrix (output size)
%  mb:                  scalar; the width of the LDPC matrix (input size)
%  d:                   string -or- scalar: degree distribution of LDPC
%  omit_hash_indices:   1*? vector; list of indices for not to hash, ie.
%                       these columns will be columns of 0,
%                       *** Must be sorted in ascending order, or will 
%                           cause MATLAB to crash (MEX-C seg fault). ***
%  ldpc_option:         struct; options for generating the LDPC
%      evencol:             bool: attempt to make the weights even
%      no4cyc:              bool: attempt to get rid of 4 cycles
%      verbost:             bool: print details of LDPC generation
%
% returns:
%  H:                   kb*mb sparse matrix; the LDPC matrix

if isnumeric(d)
    t = num2str(d);
else
    t = d;
end

H = generate_ldpc(kb, mb, t, 2, sum(100*clock), omit_hash_indices, ldpc_option.evencol, ldpc_option.no4cyc, ldpc_option.verbose);
