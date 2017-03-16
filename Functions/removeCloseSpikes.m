function [locs] = removeCloseSpikes(locs)

diffS = diff(locs) < 5;
diffS(end+1) = 0;
locs(diffS) = [];
end