function cf = crest_factor(x)

rms = sqrt(mean(x.^2));
peak = max(abs(x));
cf = peak / rms;
