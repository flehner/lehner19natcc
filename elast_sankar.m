function elasticity_sankar=elast_sankar(r,p)
% calculates precipitation elasticity using water year data
% usage: elast_sankar(r,p)
% - r = runoff time series
% - p = precipitation time series
% Time series need to be same length

xr = size(r);
xp = size(p);
if xr(1) > xr(2)
  r = r';
end
if xp(1) > xp(2)
  p = p';
end
elasticity_sankar = median( ((r-mean(r))./(p-mean(p))) * (mean(p)/mean(r)) );

return
