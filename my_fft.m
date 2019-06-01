function [f,sf] = my_fft(st,fs)
    sf = fft(st);
    sf = fftshift(sf)./(length(st)-1);
    f = [-fs/2:fs/(length(sf)-1):fs/2];
end
