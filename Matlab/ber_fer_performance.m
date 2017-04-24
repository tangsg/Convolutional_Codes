%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Viterbi Decoding using Path survivor          %
%              -----------------                  %
%           By Dama Sreekanth
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Specifications
% 1. 1/3 code rate
% 2. memory size 6
% 3. polynomials [1011011] & [1111001] [1110101]
%for ii = 0:127
%        
%    code(ii+1,1:3) = mod(sum(bitand(kron(de2bi(ii,7,'left-msb'),ones(3,1)),g),2),2).';
%    end

%%   Theory

%%
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