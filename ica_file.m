load('EEG_frontal.mat')

% ICA algorithm
addpath("C:\Users\XAVIER\Documents\MATLAB\FastICA_2.5\FastICA_25") %add package to path

[icasig] = fastica(EEG);

vec_1 = icasig(1,:);
vec_2 = icasig(2,:);
vec_3 = icasig(3,:);
vec_4 = icasig(4,:);
vec_5 = icasig(5,:);
vec_6 = icasig(6,:);
vec_7 = icasig(7,:);
vec_8 = icasig(8,:);
vec_20 = icasig(9,:);
x = [1:2193920];

% Signal using EEG signal
figure(1)
plot(x,EEG(1,:))
title('EEG  1')
figure(2)
plot(x,EEG(2,:))
title('EEG 2')
figure(3)
plot(x,EEG(3,:))
title('EEG 3')
figure(4)
plot(x,EEG(4,:))
title('EEG 4')
figure(5)
plot(x,EEG(5,:))
title('EEG 5')
figure(6)
plot(x,EEG(6,:))
title('EEG 6')
figure(7)
plot(x,EEG(7,:))
title('EEG 7')
figure(8)
plot(x,EEG(8,:))
title('EEG 8')
figure(9)
plot(x,EEG(9,:))
title('EEG 9')

%Using ICA signal
figure(1)
plot(x,vec_1)
title('Subcomponent 1 for ICA')
figure(2)
plot(x,vec_2)
title('Subcomponent 2 for ICA')
figure(3)
plot(x,vec_3)
title('Subcomponent 3 for ICA')
figure(4)
plot(x,vec_4)
title('Subcomponent 4 for ICA')
figure(5)
plot(x,vec_5)
title('Subcomponent 5 for ICA')
figure(6)
plot(x,vec_6)
title('Subcomponent 6 for ICA')
figure(7)
plot(x,vec_7)
title('Subcomponent 7 for ICA')
figure(8)
plot(x,vec_8)
title('Subcomponent 8 for ICA')
figure(9)
plot(x,vec_20)
title('Subcomponent 9 for ICA')

save('icasig.mat','icasig')