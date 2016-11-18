function [ x ] = hash_ldpc( z, H )
% computes the compressed sequence through hashing
%
% arguments:
%  z:   mb*1 vector; the translated sequence, to be hashed
%  H:   kb*mb sparse matrix; the LDPC matrix
%
% returns:
%  x:   kb*1 vector; the compressed vector

x = mod(H*z, 2);

