%% Filters
load('icasig.mat')

for i=1:9
    figure(i)
    plot(icasig(i,:))
    title(['ICA #' num2str(i)])
end

%Best ICA's: 1,2,3,5
ica_best = icasig([1,2,3,5],:);
load('ica_best.mat')
x = 1:2193920;

%% ==================LOWPASS========================%
fiir_LP = load('low_pass.mat');

G_iir_LP = fiir_LP.G;
SOS_iir_LP = fiir_LP.SOS;

[b_iir_LP,a_iir_LP] = sos2tf(SOS_iir_LP,G_iir_LP);

yiir_LP = filter(b_iir_LP,a_iir_LP,ica_best,[],2);

%yiir_LP = yiir_LP';
for i=1:4
    figure;plot(x,yiir_LP(i,:))
    xlabel('Time (s)')
    ylabel('Amplitude')
    title(['ICA #' num2str(i) 'After LowPass IIR'])   
end
for j=1:4
    [f,signal_filt] = my_fft(yiir_LP(j,:),256);
    figure
    plot(f,abs(signal_filt))
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
end

%% ==================HIGHPASS========================%
fiir_HP = load('high_pass.mat');

G_iir_HP = fiir_HP.G;
SOS_iir_HP = fiir_HP.SOS;

[b_iir_HP,a_iir_HP] = sos2tf(SOS_iir_HP,G_iir_HP);

yiir_HP = filter(b_iir_HP,a_iir_HP,yiir_LP,[],2);

%yiir_HP = yiir_HP';
for i=1:4
    figure;plot(x,(yiir_HP(i,:)))
    xlabel('Time (s)')
    ylabel('Amplitude')
    title(['ICA #' num2str(i) 'After HighPass IIR'])   
end
for j=1:4
    [f,signal_filt] = my_fft(yiir_HP(j,:),256);
    figure
    plot(f,abs(signal_filt))
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
end
save('filtered_signal.mat','yiir_HP')

%% ===============Spectrum Analyzer===============%

for i = 1:4
    [f, sf] = my_fft(ica_best(i,:),256);
    [f,sf_filtered] = my_fft(yiir_HP(i,:),256);
    figure(i)
    subplot(2,1,1)
    plot(f,abs(sf))
    title(['Signal before filtering #' num2str(i)])
    subplot(2,1,2)
    plot(f,abs(sf_filtered))
    title(['Filtered FFT Signal #' num2str(i)])
end

