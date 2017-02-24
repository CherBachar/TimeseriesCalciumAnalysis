load('spike.mat');
figure;
for i=1:length(normalSpike)
   subplot(2,2,i);
   time = length(normalSpike{i});
   plot(1:time, normalSpike{i})
   hold on;
end
hold off;
figure;
for i=1:length(fastSpike)
   subplot(2,2,i);
   time = length(fastSpike{i});
   plot(1:time, fastSpike{i})
   hold on;
end
hold off;
figure;
for i=1:length(flatSpike)
   subplot(2,2,i);
   time = length(flatSpike{i});
   plot(1:time, flatSpike{i})
   hold on;
end

%% create spike vectors
vNormalSpike = zeros(4,25);
vNormalSpike(1,1:13) = normalSpike{1}; 
vNormalSpike(1,14:end) = vNormalSpike(1,13); 
vNormalSpike(2,1:24) = normalSpike{2};
vNormalSpike(2,25) = vNormalSpike(2,24);
vNormalSpike(3,1:12) = normalSpike{3};
vNormalSpike(3,13:end) = vNormalSpike(3,12); 
vNormalSpike(4,1:24) = normalSpike{4};
vNormalSpike(4,25) = vNormalSpike(4,24);

meanNormalSpike = mean(vNormalSpike,1);
figure;plot(meanNormalSpike);
%%

vfastSpike = zeros(4,25);
vfastSpike(1,4-1:6-1) = fastSpike{1}; 
vfastSpike(1,7-1:end) = 0; 
vfastSpike(2,3-1:6-1) = fastSpike{2};
vfastSpike(2,7-1:end) = 0;
% vfastSpike(3,1:10) = fastSpike{3};
% vfastSpike(3,11:end) = 142; 
vfastSpike(4,3-1:8-1) = fastSpike{4};
vfastSpike(4,9-1:end) = 0;

meanFastSpike = mean(vfastSpike,1);
figure;plot(meanFastSpike);
%%
vFlatSpike = zeros(4,25);
vFlatSpike(1,1+1:14+1) = flatSpike{1}; 
vFlatSpike(1,15+1:end) = vFlatSpike(1,16); 
vFlatSpike(2,1+1:16+1) = flatSpike{2};
vFlatSpike(2,17+1:end) = vFlatSpike(2,18);
vFlatSpike(3,1+1:17+1) = flatSpike{3};
vFlatSpike(3,18+1:end) = vFlatSpike(3,19); 
vFlatSpike(4,1:16) = flatSpike{4};
vFlatSpike(4,17:end) = vFlatSpike(4,17);

meanFlatSpike = mean(vFlatSpike,1);
figure;plot(meanFlatSpike);

averageSpike = mean([meanNormalSpike;meanFastSpike;meanFlatSpike],1);
figure;plot(averageSpike);
%% fast spike moved to the middle

vfastSpike = zeros(4,25);
vfastSpike(1,4+7:6+7) = fastSpike{1}; 
vfastSpike(1,7+7:end) = 0; 
vfastSpike(2,3+7:6+7) = fastSpike{2};
vfastSpike(2,7+7:end) = 0;
vfastSpike(3,1+7:10+7) = fastSpike{3};
vfastSpike(3,11+7:end) = 0; 
vfastSpike(4,3+7:8+7) = fastSpike{4};
vfastSpike(4,9+7:end) = 0;

meanFastSpike = mean(vfastSpike,1);
figure;plot(meanFastSpike);


%%
save('spikeShapes.mat', 'meanNormalSpike', 'meanFastSpike', 'meanFlatSpike', 'averageSpike');