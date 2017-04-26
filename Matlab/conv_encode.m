% /*
% Convolutional encoder
%     Copyright (C) 2017  sreekanth
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  * Author: sreekanth dama
%  * Contact: sreekanth@iith.ac.in
%  **/
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

