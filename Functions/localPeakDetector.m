function [locsTemp] = localPeakDetector(locs, handles)
df_fixedF0 = handles.df_fixedF0;
n = handles.n;
localMax = 10;
locsTemp = zeros(1,length(locs));
time = handles.time;

for l = 1:length(locs)
    x1 = locs(l)-localMax;
    x2 = locs(l)+localMax;
    if (locs(l)+localMax) > time
        x2 = time;
    end
    if (locs(l)-localMax) < 1
        x1 = 1;
    end
    [~,ind] = max(df_fixedF0(n,x1:x2));
    locsTemp(l) = locs(l) - localMax - 1 + ind;
end


end