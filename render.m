close all;

figure
ax = axes;

timesteps = round(length(t)/5):100:length(t);

M(length(timesteps)) = struct('cdata',[],'colormap',[]);

plot(ax,wall(:,1),wall(:,2),'LineWidth',1)
set(ax,'XLim',xbounds,'YLim',ybounds)
grid on
set(gcf,"Position",[400 200 1100 1000])
h = animatedline(ax,'Color',[0 0 0],'LineStyle','none','Marker','.','MarkerSize',6);
for i = 1:length(timesteps)
    clearpoints(h)
    addpoints(h,r_t(:,1,timesteps(i)),r_t(:,2,timesteps(i)))
    drawnow

    M(i) = getframe(gcf);
end

vid = VideoWriter('simrender.mp4','MPEG-4');
vid.FrameRate = 60;
open(vid)
writeVideo(vid,M)
close(vid)