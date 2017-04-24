function [ y1 ] = conv_encode( ip )
%convolution encoder
%   Detailed explanation goes here
%g = [1 0 1 1 0 1 1;1 1 1 1 0 0 1;1 1 1 0 1 0 1];
% convolutional coding, rate - 1/3, generator polynomial - [133,171,165] octal
   cip1 = mod(conv(ip,[1 0 1 1 0 1 1]),2);
   cip2 = mod(conv(ip,[1 1 1 1 0 0 1]),2);
   cip3 = mod(conv(ip,[1 1 1 0 1 0 1]),2);
   cip = [cip1;cip2;cip3];
   y1 = cip(:).';

end

