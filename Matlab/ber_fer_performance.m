% /*
% Script to evaluate convolution codes
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
clear all
clc
tic
N = 10^5;  % 4-->2min 6 --> 3.5hrs
Eb_N0_dB = [3 3.5 4 4.5 5 5.5 6 6.5]-2; % multiple Eb/N0 values
%Eb_N0_dB = [5.5 6 6.5]-2; % multiple Eb/N0 values
Ec_N0_dB = Eb_N0_dB - 10*log10(3);
nErrViterbi1 = zeros(1,length(Eb_N0_dB));
nErrViterbi2 = zeros(1,length(Eb_N0_dB));

for zz=1:length(Eb_N0_dB)
zz
    ip1 = rand(1,N)>0.5;
    ip =ip1+zeros(1,N);

% normal conv encoder
[ y1 ] = conv_encode( ip );

txsym1 = 1-2*y1; %BPSK

%% channel 
n = (1/sqrt(2))*(randn(1,size(txsym1,2))+1i*randn(1,size(txsym1,2)));
txsym = txsym1 + 10^(-Ec_N0_dB(zz)/20)*n;


%% Receiver
   %% Soft Decission Decoding
   [ ipHat_v1 ] = soft_viterbi( real(txsym) );
   %% Hard decision
   [ ipHat_v2 ] = soft_viterbi( 1-2*(real(txsym)<0) );
   % counting the errors
   nErrViterbi1(zz) = size(find([ip- ipHat_v1(1:N)]),2);
   nErrViterbi2(zz) = size(find([ip- ipHat_v2(1:N)]),2);
   
end

%%
simBer_Viterbi2 = nErrViterbi2/N; % simulated ber - Viterbi decoding BER
simBer_Viterbi1 = nErrViterbi1/N; % simulated ber - Viterbi decoding BER

theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber uncoded AWGN

close all
figure
semilogy(Eb_N0_dB,theoryBer,'LineWidth',2);
 hold on
 semilogy(Eb_N0_dB,simBer_Viterbi1,'LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer_Viterbi2,'LineWidth',2);
 grid on
legend('theory - uncoded','simulation - Viterbi soft (rate-1/3,7)','simulation - Viterbi Hard (rate-1/3,7)');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for BCC with Viterbi decoding for BPSK in AWGN');
toc
%save results