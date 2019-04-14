function savefigure(figname,pos)

if nargin < 2 || isempty(pos); pos = [1,41,1920,963]; end

set(gcf,'Position',pos,'PaperPositionMode','auto');
drawnow;
pause(0.1);

saveas(gcf,[figname '.fig']);
try
    saveas(gcf,[figname '.svg']);
    saveas(gcf,[figname '.png']);
catch
    warning(['Could not save figure ''' figname '''.']);
end

end