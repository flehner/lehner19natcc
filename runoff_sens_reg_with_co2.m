function [a1,a2,b1,b2,c1,c2,d1,d2,stats]=runoff_sens_reg_with_co2(r,p,t,co2)
% calculates precipitation elasticity and temperature sensitivity using water year data
% usage: runoff_sens_reg(r,p,t)
% - r = runoff time series (% of mean, centered on 0)
% - p = precipitation time series (% of mean, centered on 0)
% - t = temperature time series
% - co2 = CO2 anomaly to 300ppm
% Time series need to be same length
% Returns (in this order):
% - (1) P sensitivity
% - (2) P sensitivity uncertainty
% - (3) T sensitivity
% - (4) T sensitivity uncertainty
% - (5) CO2 sensitivity
% - (6) CO2 sensitivity uncertainty
% - (7) PT sensitivity
% - (8) PT sensitivity uncertainty

% -- regression parameter
add_one = 0; % 0=no, 1=yes

xr = size(r);
xp = size(p);
xt = size(t);
xc = size(co2);
if xr(1) < xr(2)
  r = r';
end
if xp(1) < xp(2)
  p = p';
end
if xt(1) < xt(2)
  t = t';
end
if xc(1) < xc(2)
  co2 = co2';
end
if add_one == 0
  X   = [p t co2 (p.*t)];
else
  X   = [ones(size(p)) p t co2 (p.*t)];
end
[tmp ci] = regress(r,X);
a1 = tmp(1+add_one);
b1 = tmp(2+add_one);
a2 = ci(1+add_one,:);
b2 = ci(2+add_one,:);
c1 = tmp(3+add_one);
c2 = ci(3+add_one);
d1 = tmp(4+add_one); % interaction term
d2 = ci(4+add_one); % CI on interaction term


return
