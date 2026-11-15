function v = p95(x)
%P95  95th percentile (toolbox-free).
x = sort(x(~isnan(x)));
if isempty(x), v = NaN; return, end
v = x(max(1, ceil(0.95*numel(x))));
end
