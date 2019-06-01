% Data loading
load('wavelet_dec_final.mat')

% Creating a 2D dimension with zeros
X = zeros(size(wavelet_dec_final,2),20);

% Labeling
load('G8Data\zeus_data_pat5_v2.mat')
BIS = data.BIS;

label = ones(4283,1);  % 1 - awake
label(1800:3990) = 2;  % 2 - asleep

% Transforming wavelet_ced_final into a 2D matriz
X(:,1:5)   = wavelet_dec_final(1,:,:);
X(:,6:10)  = wavelet_dec_final(2,:,:);
X(:,11:15) = wavelet_dec_final(3,:,:);
X(:,16:20) = wavelet_dec_final(4,:,:);

% PCA
[coeff,score,latent,tsquared,explained] = pca(X,'NumComponents',3);
[coeff2,score2,latent2,tsquared2,explained2] = pca(X,'NumComponents',2);

%plot(score2)
%plot(score)

save('score_pc_3d.mat','score')
save('score_pca_2d.mat','score2')

figure(1)
scatter3(score(label==1,1),score(label==1,2),score(label==1,3),'g')
hold on 
scatter3(score(label==2,1),score(label==2,2),score(label==2,3),'r')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Green = awake', 'Red = asleep')
title('PCA 3D')

figure(2)
scatter(score2(label==1,1),score2(label==1,2),'g')
hold on 
scatter(score2(label==2,1),score2(label==2,2),'r')
legend('Green = awake', 'Red = asleep')
xlabel('x')
ylabel('y')
title('PCA 2D')

% Dissimilarity matrix
D = pdist(X,'euclidean'); 
% MDS
[Y,eigvals] = cmdscale(D);
% save('mds.mat','Y')

figure(3)
scatter(Y(label==1,1),Y(label==1,2),'g')
hold on
scatter(Y(label==2,1),Y(label==2,2),'r')
legend('Green = awake', 'Red = asleep')
xlabel('x')
ylabel('y')
title('MDS 2D')


% Dissimilarity matrix of the output
D2 = pdist(Y,'euclidean');
plot(D,D2,'bo',[0 50],[0 50],'k--')
xlabel('Dissimilarities')
ylabel('Disparities')
title('MDS')

% Another MDS 

% Distance
dist_eu=pdist(X,'euclidean');
dist_c=pdist(X,'cityblock');
dist_cheb=pdist(X,'chebychev');

% Matrix
matrix_disteu=squareform(dist_eu);
matrix_distc=squareform(dist_c);
matrix_distcheb=squareform(dist_cheb);

% Using chevchev Distance - 2D
[y2D_c, stress2D_c, disparities2D_c] = mdscale(matrix_disteu, 2, 'criterion', 'metricsstress');
%[y2D_c_1, stress2D_c_1, disparitiesD_c_1] = mdscale(matrix_distcheb, 2, 'criterion', 'stress');

% Using chevchev Distance - 3D
[y3D_che, stress3D_che, disparities3D_c] = mdscale(matrix_disteu, 3, 'criterion', 'metricsstress');
%[y3D_che_1, stress3D_che_1, disparities3D_che_1] = mdscale(matrix_distcheb, 3, 'criterion', 'stress');

y_data = y3D_che(:,1);
y_data_2 = y3D_che(:,2);

%load y_data;
%load y_data_2;

plot(y_data)
hold on 
plot(y_data_2)