function [locsTemp] = findLocalMax(locs, df_fixedF0, n, time, localMax)
locsTemp = zeros(1,length(locs));
time = size(df_fixedF0,2);

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
    
    %make sure it's between 1-time
    if locsTemp(l) < 1
        locsTemp(l) = 1;
    end
    if locsTemp(l) > time
        locsTemp(l) = time;
    end        
end

[locsTemp] = checkLocalMax(locsTemp, df_fixedF0, n, time);
[locsTemp] = removeCloseSpikes(locsTemp);
end

function [locs] = checkLocalMax(locs, df_fixedF0, n, time)
toRemove = [];
i=1;
for l = 1:length(locs)
    x1 = locs(l)-1;
    x2 = locs(l)+1;
    if x2 > time
        x2 = time;
    end
    if x1 < 1
        x1 = 1;
    end
   if  (df_fixedF0(n,locs(l)) <  df_fixedF0(n,x2)) ...
           || (df_fixedF0(n,locs(l)) <  df_fixedF0(n,x1))
        toRemove(i) = l;
        i = i+1;
   end
end
if ~isempty(toRemove)
    locs(toRemove) = [];
end
end

