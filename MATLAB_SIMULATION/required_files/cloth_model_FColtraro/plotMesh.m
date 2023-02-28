function [super,curve] = plotMesh(X,T,nodos_borde)
figure(1)
super = trisurf(T,X(:,1),X(:,2),X(:,3),'edgecolor', 'k','facecolor','w');
hold on
gamma0 = X(nodos_borde,:);
curve = plot3(gamma0(:,1),gamma0(:,2),gamma0(:,3),'.','color','k','MarkerSize',10);
light               % add a light
lighting gouraud    % preferred lighting for a curved surface
camzoom(2)          % zoom into scene
axis off
axis equal
axis([-0.4, 0.4, 0, 1, -0.1, 0.4])
%axis tight
% set axis equal and remove axis
xlabel('X'); ylabel('Y');zlabel('Z');
end