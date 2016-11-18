function [H] = ldpc_generate_n4c(m,n,d,~,seed,dope_indices,evencol,no4cycle,verbose)
% MEX file: generates sparse LDPC parity check matrix over GF2
% 
% arguments:
%  m:               scalar; number of parity checks
%  n:               scalar; blocklength
%  d:               string: degree distribution
%  ~                GF base (only 2 now)
%  seed:            number; initializes random generator
%  dope_indices:    1*? vector; indices not to hash
%  evencol:         bool: 1 for only col, 0 for both col and row
%  no4cycle:        bool: 1 for no4cycle, 0 to allow 4 cycles
%  verbose:         bool: 1 for verbose
%
% returns:
%  H:               m*n sparse matrix; the LDPC matrix
