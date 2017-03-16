function [handles] = findCaSynchrony(trace, handles)
%% Initialize variables

cells = trace';
time_points = size(cells,1);
num_cells = size(cells,2);
%time = 1:1:time_points;

window = [];
noverlap = [];
nfft = [];
fs = 0.3; %frames per second

pNetAv = [];
pRandAv = [];
%% Cell power spectral density
%calculate the cell power spectral density


for i = 1:num_cells
    [Pi, f ]= pwelch(cells(:,i), window, noverlap, nfft, fs );

    if i == 1
        Pcell = Pi;
        f1 = f;
    else
        Pcell = Pcell + Pi;
        f1 = f1 + f;
    end    
end

Pcell = Pcell / num_cells;
f1 = f1 / num_cells;

%% Network power spectral density
%calculate the network power spectral density
for n = 1:time_points
    Pmean = mean(cells,2);
    [Pnet, freq] = pwelch(Pmean, window, noverlap, nfft, fs);

    %% Randomize input
    %Iterate over number of time points to accumulate power and phase randomize 
    
    %randIndex = randperm(size(cells,1)*size(cells,2));
    %randCells = reshape(cells(randIndex),size(cells));
    %Prandmean = mean(randCells,2);
    
    PrandNet = phaseran(Pmean, num_cells);%phase random over num_cells blocks
    PrandNet = mean(PrandNet,3); %get mean for all cells
    [Prand, frand] = pwelch(PrandNet, window, noverlap, nfft, fs);
    
    %Accumulate network and phaserand power
    if isempty(pNetAv)
        pNetAv = Pnet; 
        pRandAv = Prand;
    else
        pNetAv = Pnet + pNetAv; 
        pRandAv = Prand + pRandAv;
    end
    
    synchTemp(:,n) = (Pnet-Prand)./(Pcell-Prand);
end

pNetAv = pNetAv / (time_points);
pRandAv = pRandAv / (time_points);

%Calculate synchronicity from the network and randphase average
synch1 = (pNetAv-pRandAv)./(Pcell-pRandAv);
%% Plot 

%Spectral analysis
figure, loglog(freq,Pcell,'r');  
title('Spectral analysis of Ca2+ time series data', 'FontSize',12.5);
xlabel('Frequency (Hz)');
ylabel('[\DeltaF/F]^2/Hz');
hold on;
loglog(freq, pNetAv,'b');
loglog(freq, pRandAv,'c');
legend('power cell', 'power network', 'random');
hold off;

%Synchronicity
figure;
synch2 = mean(synchTemp, 2);
semilogx(freq, synch1);
title('Synchronicity', 'FontSize',12.5);
xlabel('Frequency (Hz)');
ylabel('[Network - PhaseRand]/[Cell - PhaseRand] Power Ratio (= Synchronicity)');
%% Save variables

handles.synch.synch = synch1;
handles.synch.freq = freq;
handles.synch.Pcell = Pcell;
handles.synch.Pnet = Pnet;
handles.synch.Prand = Prand;

end