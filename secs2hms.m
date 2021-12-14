function hms = secs2hms(secs)
% Converts a number of secs to a string in hh:mm:ss format

hms = sprintf('%02d:%02d:%02d', ...
    floor(secs/3600), ... % hours
    floor(mod(secs/60,60)), ... % minutes
    floor(mod(secs,60)) ... % seconds
    );

end