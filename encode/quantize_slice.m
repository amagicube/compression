function [ u ] = quantize_slice( s, Q, Q0, b )
% quantizes the source
%
% arguments:
%  s:   n*1 vector; sample to be quantized
%  Q:   m*n sparse matrix; the quantization matrix
%  Q0:  m*1 vector; the offset vector
%  b:   scalar; the size of the translator output
%
% returns:
%  u:   m*1 vector; the quantized sequence

num_bins = 2^b;

u = min(max(floor(Q*s + Q0),-num_bins/2),num_bins/2-1);
