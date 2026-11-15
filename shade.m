function shade(ax, S, yl)
%SHADE  Grey band over the radar measurement-loss intervals of scenario S.
for k = 1:size(S.blackout, 1)
    t0 = S.blackout(k, 1);  t1 = t0 + S.blackout(k, 2);
    patch(ax, [t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], [0.88 0.88 0.88], 'EdgeColor', 'none');
end
end
