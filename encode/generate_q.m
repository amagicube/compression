function [ Q, Q0 ] = generate_q( n, m, q_option )
% generate the quantization matrix
%
% arguments:
%  n:           scalar; width of the quantization matrix
%  m:           scalar; height of the quantization matrix
%  q_option:    struct; contains information about generating Q
%
% returns:
%  Q:           m*n sparse matrix; the quantization matrix
%  Q0:          m*1 vector; the offset vector

Q = speye(n)/q_option.w;
Q0 = zeros(m,1);