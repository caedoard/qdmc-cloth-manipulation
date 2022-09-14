function hacerPeli(phiPositions,T,nodos_borde,peli,vel_vid) %,x,y,z,x2,y2,z2)
%triangulamos para hacer el plot
[Xtri,Ttri] =  traingulateQuadMesh(phiPositions{1},T); 
[super,curve] = plotMesh(Xtri,Ttri,nodos_borde);

n_frames = size(phiPositions,2);
clear VID
if ~isempty(peli)
    set(gca, 'nextplot', 'replacechildren');
    VID(1) = getframe;
end    
%actualizamos la figura
for i=2:vel_vid:n_frames
    Xnew = phiPositions{i};
    [super,curve] = updatePlot(super,curve,Xnew,T,nodos_borde);
%     hold on
%     plot3(x(i),y(i),z(i),'.r')
%     hold on
%     plot3(x2(i),y2(i),z2(i),'.r')
    pause(0.01)
%     w = waitforbuttonpress;
    if ~isempty(peli)
        if vel_vid == 1
           VID(i) = getframe;
        else   
           VID(floor(i/vel_vid)+1) = getframe;
        end
    end 
end 

if ~isempty(peli)
    %guardamos el video
    v = VideoWriter(peli);
    v.Quality = 95;
    v.FrameRate = 100;
    open(v)
    writeVideo(v,VID)
    close(v)
end    