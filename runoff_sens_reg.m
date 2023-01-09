function [a1,a2,b1,b2,c1,c2,stats]=runoff_sens_reg(r,p,t)
% calculates precipitation elasticity and temperature sensitivity using water year data
% usage: runoff_sens_reg(r,p,t)
% - r = runoff time series (% of mean, centered on 0)
% - p = precipitation time series (% of mean, centered on 0)
% - t = temperature time series
% Time series need to be same length
% Returns (in this order):
% - (1) P sensitivity
% - (2) P sensitivity uncertainty
% - (3) T sensitivity
% - (4) T sensitivity uncertainty

% -- regression parameter
add_one = 0; % 0=no, 1=yes

xr = size(r);
xp = size(p);
xt = size(t);
if xr(1) < xr(2)
  r = r';
end
if xp(1) < xp(2)
  p = p';
end
if xt(1) < xt(2)
  t = t';
end
if add_one == 0
  X   = [p t (p.*t)];
else
  X   = [ones(size(p)) p t (p.*t)];
end
[tmp ci] = regress(r,X);
a1 = tmp(1+add_one);
b1 = tmp(2+add_one);
a2 = ci(1+add_one,:);
b2 = ci(2+add_one,:);
c1 = tmp(3+add_one); % interaction term
c2 = ci(3+add_one); % CI on interaction term


return
