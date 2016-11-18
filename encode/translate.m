function [ z ] = translate( u, b, bin_code )
% translates the quantized sequence into a 0-1 sequence
%
% arguments:
%  u:           m*1 vector; the vector to be translated
%  b:           scalar: the size of the translator output
%  bin_code:    string; the translator choice: 'bin' for standard binary
%               code, 'pam' for PAM gray code, 'qam' for QAM gray code
%
% returns:
%  z:           mb*1 vector; the translated sequence

num_bins = 2^b;
u_adj = u + num_bins/2;
if strcmp(bin_code, 'bin')
    z = reshape(de2bi(u_adj, b,'left-msb')', [],1);
else
    z = reshape(de2bi(bin2gray(u_adj,bin_code,num_bins), b,'left-msb')', [],1);
end