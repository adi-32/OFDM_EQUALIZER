clc;
clear all;
SNR_dB=0:2:24;
SNR=10.^(SNR_dB/10);
N=128;
N_iterations=1e3;
FFT_length=N;
CP_length=16;

Mod_type='BPSK'; %BPSK, QPSK, 16QAM, 64QAM, or 256QAM
Equalizer_type='LMMSE';%LZF, or LMMSE

switch Mod_type
    case {'BPSK'}
        Nb=1;
        BPSK_Mod=comm.PSKModulator(2, 'BitInput', true, 'PhaseOffset', 0);
        BPSK_Demod=comm.PSKDemodulator(2, 'BitOutput', true, 'PhaseOffset', 0);        
    case {'QPSK'}
        Nb=2;
        QPSK_Mod=comm.PSKModulator(4, 'BitInput', true, 'PhaseOffset', pi/4);
        QPSK_Demod=comm.PSKDemodulator(4, 'BitOutput', true, 'PhaseOffset', pi/4);
    case {'16QAM'}
        Nb=4;
        QAM16_Mod=comm.RectangularQAMModulator(16, 'BitInput', true,....
            'NormalizationMethod','Average power','SymbolMapping','Gray');
        QAM16_Demod=comm.RectangularQAMDemodulator(16, 'BitOutput', true,....
            'NormalizationMethod','Average power','SymbolMapping','Gray');
    case {'64QAM'}
        Nb=6;
        QAM64_Mod=comm.RectangularQAMModulator(64, 'BitInput', true,....
            'NormalizationMethod','Average power','SymbolMapping','Gray');
        QAM64_Demod=comm.RectangularQAMDemodulator(64, 'BitOutput', true,....
            'NormalizationMethod','Average power','SymbolMapping','Gray');        
    case {'256QAM'}
        Nb=8;
        QAM256_Mod=comm.RectangularQAMModulator(256, 'BitInput', true,....
            'NormalizationMethod','Average power','SymbolMapping','Gray');
        QAM256_Demod=comm.RectangularQAMDemodulator(256, 'BitOutput', true,....
            'NormalizationMethod','Average power','SymbolMapping','Gray');
end


for n=1:length(SNR)
    
for k=1:N_iterations
    h=jake(120, 6);
    H=diag(fft(h, N));
    
    switch Equalizer_type
        case {'LZF'}
          C=inv(H'*H)*H';              
        case {'LMMSE'}
          C=inv(H'*H+(1/SNR(n)*eye(N)))*H';  
    end
d=round(rand(N*Nb, 1));    

switch Mod_type
    case {'BPSK'}
    X=step(BPSK_Mod, d);    
    case {'QPSK'}
    X=step(QPSK_Mod, d);        
    case {'16QAM'}
    X=step(QAM16_Mod, d);        
    case {'64QAM'}
    X=step(QAM64_Mod, d);                
     case {'256QAM'}
    X=step(QAM256_Mod, d);        
end    

x=sqrt(FFT_length)*ifft(X.');  %%% bank of Modulators
x_cp=[x(N-CP_length+1:N) x];  %%% CP insertion
y=filter(h,  1,  x_cp);                        %% Rayleigh channel

y=awgn(y, SNR_dB(n),'measured'); % Add AWGN

y1=y(CP_length+1:N+CP_length); %% CP removal

Y1=1/sqrt(FFT_length)*fft(y1.');
Y=C*Y1;

switch Mod_type
    case {'BPSK'}
    Y_est=step(BPSK_Demod, Y);    
    case {'QPSK'}
    Y_est=step(QPSK_Demod, Y);        
    case {'16QAM'}
    Y_est=step(QAM16_Demod, Y);        
    case {'64QAM'}
    Y_est=step(QAM64_Demod, Y);                
     case {'256QAM'}
    Y_est=step(QAM256_Demod, Y);        
end    

error(k)=sum(Y_est~=d);

end
BER(n)=mean(error)/(N*Nb);
end
%figure;
semilogy(SNR_dB, BER,'-br','linewidth', 1.5); hold on
xlabel('SNR [dB]'); ylabel('BER');
title('SISO-OFDM, Rayleigh, N=64, N_{CP}=16')
axis([0 24 1e-4 1])
grid



