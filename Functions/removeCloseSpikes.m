function [locs] = removeCloseSpikes(locs)

diffS = diff(locs) < 3;
diffS(end+1) = 0;
locs(diffS) = [];
end