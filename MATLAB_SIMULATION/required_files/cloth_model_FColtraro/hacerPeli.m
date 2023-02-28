function hacerPeli(phiPositions,T,nodos_borde,peli,vel_vid,dt)
%triangulamos para hacer el plot
[Xtri,Ttri] =  triangulateQuadMesh(full(phiPositions{1}),T); 
[super,curve] = plotMesh(Xtri,Ttri,nodos_borde);

n_frames = size(phiPositions,2);
clear VID
if ~isempty(peli)
    set(gca, 'nextplot', 'replacechildren');
    VID(1) = getframe(gcf);
end    
%actualizamos la figura
for i=2:vel_vid:n_frames
    Xnew = triangulateQuadMesh(full(phiPositions{i}),T);
    [super,curve] = updatePlot(super,curve,Xnew,nodos_borde);
    pause(dt)
    %w = waitforbuttonpress;
    if ~isempty(peli)
        if vel_vid == 1
           VID(i) = getframe(gcf);
        else   
           VID(floor(i/vel_vid)+1) = getframe(gcf);
        end
    end 
end 

if ~isempty(peli)
    %guardamos el video
    v = VideoWriter(peli);
    v.Quality = 50;
    v.FrameRate = 100;
    open(v)
    writeVideo(v,VID)
    close(v)
end    