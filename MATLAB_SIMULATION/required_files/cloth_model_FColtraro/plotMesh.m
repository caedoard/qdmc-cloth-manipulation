function [super,curve] = plotMesh(X,T,nodos_borde)
% figure(1)
super = trisurf(T,X(:,1),X(:,2),X(:,3),'edgecolor', 'k','facecolor','w');
hold on
%gamma0 = X(nodos_borde.nodes_bnd,:);
curve = [];% plot3(gamma0(:,1),gamma0(:,2),gamma0(:,3),'.','color','k','MarkerSize',10);
% light               % add a light
% lighting gouraud    % preferred lighting for a curved surface
ylim([0,2.5]); xlim([-0.5,0.5]); zlim([-0.5,0.5])
view([-50 20])
camzoom(2.0)   % zoom into scene
% axis off
% axis tight
% axis equal
% set axis equal and remove axis
% xlabel('X'); ylabel('Y');zlabel('Z');