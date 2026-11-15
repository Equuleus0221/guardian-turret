function save_fig(fig, name, width_in, height_in)
%SAVE_FIG  Save a figure as vector PDF (and PNG) in ./figures at a fixed size.
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end
set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1), pos(2), width_in, height_in], 'Color', 'w', ...
    'PaperUnits', 'inches', 'PaperSize', [width_in, height_in], ...
    'PaperPosition', [0, 0, width_in, height_in]);
drawnow;
pdf = fullfile(outdir, [name, '.pdf']);
png = fullfile(outdir, [name, '.png']);
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, pdf, 'ContentType', 'vector');
    exportgraphics(fig, png, 'Resolution', 200);
else
    print(fig, pdf, '-dpdf');
    print(fig, png, '-dpng', '-r200');
end
end
