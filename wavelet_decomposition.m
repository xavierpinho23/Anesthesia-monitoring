%% Wavelet Decomposition

load('filtered_signal.mat')
% Separate signal in time intervals
n = size(yiir_HP,2)/256;

mother_wav = 'db4';
dec_level = 9;
wavelet_dec = [];

for i=1:4
    for j=1:2:n-5
        p_segment = yiir_HP(i,(j*256):((j+5)*256));
        wavelet_dec(i,(j+1)/2,:)=wav_coef(p_segment,mother_wav, dec_level);
    end
end

save('wavelet_dec.mat','wavelet_dec')

%% Bands

load('wavelet_dec.mat')
% delta 8,7,6
% theta 5
% alpha 4
% beta 3
% gamma 2,1

% Aply the mean to each window by band

wavelet_dec_final = [];

for i=1:4 %leeds
    for j = 1:size(wavelet_dec,2)
    wavelet_dec_final(i,j,1) = wavelet_dec(i,j,3);
    wavelet_dec_final(i,j,2) = wavelet_dec(i,j,4);
    wavelet_dec_final(i,j,3) = wavelet_dec(i,j,5);
    wavelet_dec_final(i,j,4) = wavelet_dec(i,j,6);
    wavelet_dec_final(i,j,5) = (wavelet_dec(i,j,9) + wavelet_dec(i,j,7) + wavelet_dec(i,j,8)/3);

    end
end


for i=1:5
    for j=1:4
       wavelet_dec_final(j,:,i) = hampel(wavelet_dec_final(j,:,i),100,2);
    end
end


save('wavelet_dec_final.mat','wavelet_dec_final')

titles = {'Gamma','Beta','Alpha','Theta','Delta'};

for j=1:5
    figure(j)
    suptitle(titles{j})
    subplot(2,2,1);        
    stem(wavelet_dec_final(1,:,j),'Marker','none');
    title('component 1')
    subplot(2,2,2);
    stem(wavelet_dec_final(2,:,j),'Marker','none');
    title('component 2')
    subplot(2,2,3);        
    stem(wavelet_dec_final(3,:,j),'Marker','none');
    title('component 3')
    subplot(2,2,4);
    stem(wavelet_dec_final(4,:,j),'Marker','none');
    title('component 4')
end
