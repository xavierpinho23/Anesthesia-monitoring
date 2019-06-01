% Data loading
load('score_pc_3d.mat') %PCA 3d
load('score_pca_2d.mat') %PCA 2d
load('mds.mat') %MDS

% Labeling
label = ones(4283,1);  % 1 - awake
label(1800:3990) = 2;  % 2 - asleep

% Clustering

% K-means 3D
[idx_3d,c_3d] = kmeans(score,2);
crosstab(label,idx_3d)
% K-means 2D
[idx_2d,c_2d] = kmeans(score2,2);
crosstab(label,idx_2d)
% Kmeans MDS
[idx_mds,c_mds] = kmeans(Y,2);
crosstab(label,idx_mds)

% plot clustering PCA 3D
figure(1)
plot(score(idx_3d==1,1),score(idx_3d==1,2),'r.','MarkerSize',12)
hold on
plot(score(idx_3d==2,1),score(idx_3d==2,2),'b.','MarkerSize',12)
hold on
plot(c_3d(:,1),c_3d(:,2),'ko','MarkerSize',12,'LineWidth',2);

% Plot clustering PCA 2D
figure(2)
plot(score2(idx_2d==1,1),score2(idx_2d==1,2),'r.','MarkerSize',12)
hold on
plot(score2(idx_2d==2,1),score2(idx_2d==2,2),'b.','MarkerSize',12)
hold on
plot(c_2d(:,1),c_2d(:,2),'ko','MarkerSize',12,'LineWidth',2);

% plot clustering MDS
figure(3)
plot(Y(idx_mds==1,1),Y(idx_mds==1,2),'r.','MarkerSize',12)
hold on
plot(Y(idx_mds==2,1),Y(idx_mds==2,2),'b.','MarkerSize',12)
hold on
plot(c_mds(:,1),c_mds(:,2),'ko','MarkerSize',12,'LineWidth',2);

% Hirearchical clustering - PCA 3D

cheb = pdist(score,'chebychev');
clustTreeCheb = linkage(cheb,'ward');
dendrogram(clustTreeCheb)

cophenet(clustTreeCheb,cheb)
c1 = cluster(clustTreeCheb,'maxclust',2,'Depth',1);
scatter3(score(c1==1,1),score(c1==1,2),score(c1==1,3))
hold on
scatter3(score(c1==2,1),score(c1==2,2),score(c1==2,3))
title('Hierarchical Clustering - PCA 3D')

crosstab(label,c1)

% Hirearchical clustering - PCA 2D

euc = pdist(score2,'euclidean');
clustTreeEuc = linkage(euc,'ward');
dendrogram(clustTreeEuc)

cophenet(clustTreeEuc,euc)
c2 = cluster(clustTreeEuc,'maxclust',2,'Depth',1);
gscatter(score(:,1),score(:,2),c2)
title('Hierarchical Clustering - PCA 2D')

crosstab(label,c2)

%%% hierarchical clustering MDS
euc_pca_2d = pdist(Y,'euclidean');
clustTreeEuc_pca_2d = linkage(euc_pca_2d,'ward');
dendrogram(clustTreeEuc_pca_2d)

cophenet(clustTreeEuc_pca_2d,euc_pca_2d)
c4 = cluster(clustTreeEuc_pca_2d,'maxclust',2,'Depth',1);
gscatter(Y(:,1),Y(:,2),c4)
title('Hierarchical Clustering - MDS')

crosstab(label,c4)

