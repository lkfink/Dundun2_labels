% LF new attempt at acoustic feats for dundun
% 20210801
% 20221210

% Set directories
cwd = '/Users/lauren.fink/Documents/Projects/dundun/Code/acousticAnalyses/Madita_stims';
stim_dir = fullfile(cwd, 'stims/'); %folder of stimuli 
audiofiles = dir(fullfile(stim_dir,'*.wav')); % create directory of wave files
nfiles = length(audiofiles); % how many files we have
figpath = fullfile(cwd, '/figures'); % where to save optional figures

% add path to mir toolbox
addpath(genpath('/Users/lauren.fink/Documents/MATLAB_utils/MIRtoolbox1.7.2'));
addpath(genpath('/Users/lauren.fink/Documents/MATLAB_utils/helperFunctions'));

% flag if want to use KF onsets
useKFonsets = 0;
essentia = 0;

%-------------------------------------------------------------%
% set constants and initialize variables
%-------------------------------------------------------------%
% previous paper used 10 ms time windows and 1 ms steps
winLength = 480*5; %round(0.05*fs); 2400 % 50 ms
overlapLength = winLength*(3/4); % 75% overlap
as = table; % table to store audio data
nr = 1; % row counter for our table
toplot = 0; % flag whether to plot data or not

% set RQA parameters
del = 10;
emb = 2;
rad = 0.3;
%-------------------------------------------------------------%

% loop through all audio files
% extract features we want and save to table
for ifile = 1:nfiles

    thisWavName = fullfile(stim_dir, audiofiles(ifile).name);

    % ONLY do for KF files
    % 26.11.21
    if useKFonsets
        f = sum(strcmp(audiofiles(ifile).name(1:end-4), kfstims)) > 0;
    end
%     if f <= 0 % stupid code. fix. or uncomment for kf version
%         continue
%     end

    
    [audio,fs] = audioread(thisWavName);
    audio = audio(:,1); % only take single channel
    as.stim{nr,1} = audiofiles(ifile).name(1:end-4);
    as.audio{nr,1} = audio; % save raw audio to table
    as.fs(nr,1) = fs;
    
    %-------------------------------------------------------------%
    % amplitude envelope
    [upper, lower] = envelope(audio, winLength, 'rms');
    [b,a] = butter(3, 50/(fs/2), 'low'); % filter
    filtenv = filtfilt(b,a,upper(:,1)); % filter amp env?
    normenv = zscore(filtenv);
    
    % find peaks in the amplitude envelope to define note onsets and durations
    [peaks, noteOnsets, noteDurs] = findpeaks(normenv, ...
        'MinPeakDistance', 48*40, ... % 40 ms between peaks
        'MinPeakProminence',.1, ...
        'WidthReference', 'halfprom');
    
    
    %-------------------------------------------------------------%
    % pitch
    [f0,idx] = pitch(audio,fs, ...
        'Method','PEF', ...
        'Range', [40, 600], ...
        'WindowLength',winLength/5, ... % 50 ms window. 10 ms wins
        'OverlapLength',overlapLength/5, ... % 75% overlap
        'MedianFilterLength',3); % 1 = no smoothing
    as.pitch{nr,1} = f0;
    
    
    %-------------------------------------------------------------%
    % entropy
    [we, te] = pentropy(audio, fs, 'Instantaneous', true, 'Scaled', true);
    entTime = te*fs; % convert entropy time to samples
    as.ent{nr,1} = we;
    
    %-------------------------------------------------------------%
    % new 25.11.22
    % spectral flux
    flux = spectralFlux(audio,fs, 'Window', rectwin(round(fs*0.05)), 'OverlapLength',1800); 
    as.flux{nr,1} = flux;
    as.fluxMean(nr,1) = mean(flux);

    %-------------------------------------------------------------%
    % recurrance quantification analysis
    [RP, RESULTS, PARAMETERS, b] = MDRQA(flux,emb,del,'euc',rad,1);
    as.rec_results{nr,1} = RESULTS;
    as.recurr(nr,1) = RESULTS(2);
    as.determ(nr,1) = RESULTS(3);
    as.meanL(nr,1) = RESULTS(4);
    as.maxL(nr,1) = RESULTS(5);
    as.entrL(nr,1) = RESULTS(6);

%     figure()
%     imagesc(RP*-1)
%     colormap(gray);
%     xlabel('Flux')
%     ylabel('Flux')

    % rec, det, lam probably most interesting in this context 

    %  RESULTS is a double-variable holding the following recurrence variables:
    %    1.  Size of the RP
    %    2.  %REC  - percentage of recurrent points
    %    3.  %DET  - percentage of diagonally adjacent recurrent points
    %    4.  MeanL - average length of adjacent recurrent points
    %    5.  MaxL  - maximum length of diagonally adjacent recurrent points
    %    6.  EntrL - Shannon entropy of distribution of diagonal lines
    %    7.  %LAM  - percentage of vertically adjacent recurrent points
    %    8.  MeanV - average length of diagonally adjacent recurrent points
    %    9.  MaxV  - maximum length of vertically adjacent recurrent points
    %    10. EntrV - Shannon entropy of distribution of vertical lines
    
        
    %-------------------------------------------------------------%
    % Plot our measures
    t = (0:length(audio)-1)/fs;
    x_peaks = t(noteOnsets); % get time how we need it
    tf0 = (idx - 1)/fs;
    pitchTime = tf0*fs; % convert pitch times to samples
    
    if toplot
        figure()
        subplot(3,1,1)
        plot(t,normenv,  x_peaks,peaks,'pg')
        hold on
        plot(t, audio)
        xlabel('Time (s)')
        ylabel('Amplitude')
        legend('Normalized Envelope','Note Onsets', 'Original Signal')
        title('Amplitude Envelope & Note Dectection')
        hold on
        
        subplot(3,1,2)
        
        yyaxis left; plot(t,audio)
        ylabel('Amplitude')
        yyaxis right; plot(tf0,f0)
        xlabel('Time (s)')
        ylabel('Pitch (Hz)')
        ylim([40 600])
        title('Pitch extraction')
        hold on
        
        % plot entropy
        subplot(3,1,3)
        plot(t,audio)
        hold on
        plot(te, we)
        ylabel('Spectral Entropy')
        xlabel('Time (sec)')
        title('Wiener Entropy')
        suptitle(audiofiles(ifile).name(1:end-4))
    end
    
    %-------------------------------------------------------------%
    % Get features by note (rather than continuous)
    % take note length as note onset +20 ms
    notedursamps = fs/(1000/20);
    noteOffsets = noteOnsets + notedursamps;

    % Here NEW 26.11.21
    % If this is a stimulus we have onsets from Klaus, use klaus onsets
    if useKFonsets
        if sum(strcmp(as.stim{nr,1}, kf.stim) > 0) % we have data from KF
            if essentia
                nm = strcmp(kf.type, 'essentia');
            else
                nm = strcmp(kf.type, 'manual');
            end
            kfmask = strcmp(as.stim{nr,1}, kf.stim);
            compmask = kfmask & nm;
            noteOnsets = kf.onset(compmask);
            noteOnsets = noteOnsets*fs;
            noteOffsets = noteOnsets + notedursamps;
            %badinds =  find(noteOffsets(end) > length(we)

        end
    end
    % Just redo this here ^ and save separately to compare to formerly
    % calculated table.

    for i = 1:length(noteOnsets)
        
        % amplitude
        amp_means(i) = nanmean(normenv(noteOnsets(i):noteOffsets(i)));
        
        % pitch
        pitchMask = pitchTime >= noteOnsets(i) & pitchTime <= noteOffsets(i);
        pitchVals = f0(pitchMask);
        pitch_means(i) = nanmean(pitchVals);
        
        % entropy
        entMask = entTime >= noteOnsets(i) & entTime <= noteOffsets(i);
        ent_means(i) = nanmean(we(entMask));

        % todo flux

    end
    
    % save values and differences for each note
    % amp
    as.env_notes{nr,1} = amp_means;
    as.env_noteMean{nr,1} = nanmean(amp_means);
    as.env_noteDiffs{nr,1} = diff(amp_means);
    as.env_noteMeanDiff{nr,1} = nanmean(abs(as.env_noteDiffs{nr,1}));
    
    % pitch
    as.pitch_notes{nr,1} = pitch_means;
    as.pitch_noteMean{nr,1} = nanmean(pitch_means);
    as.pitch_noteDiffs{nr,1} = diff(pitch_means);
    as.pitch_noteMeanDiff{nr,1} = nanmean(abs(as.pitch_noteDiffs{nr,1}));
    
    % ent
    as.ent_notes{nr,1} = ent_means;
    as.ent_noteMean{nr,1} = nanmean(ent_means);
    as.ent_noteDiffs{nr,1} = diff(ent_means);
    as.ent_noteMeanDiff{nr,1} = nanmean(abs(as.ent_noteDiffs{nr,1}));
    
    %-------------------------------------------------------------%
    % Calculate IOIs
    % note onsets are in samples. convert to ms
    notems = noteOnsets/fs;
    as.IOI_notes{nr,1} = diff(notems)*1000;
    as.IOI_noteMean(nr,1) = nanmean(as.IOI_notes{nr,1});
    as.IOI_noteDiffs{nr,1} = diff(as.IOI_notes{nr,1});
    as.IOI_noteMeanDiff(nr,1) = nanmean(abs(as.IOI_noteDiffs{nr,1}));
    
    % Ratio
    % interval1/(interval1+interval2)
    IOIs = as.IOI_notes{nr,1};
    for i = 1:length(IOIs) -1 % can't do last because don't have following interval
        ratio(i) = IOIs(i) / (IOIs(i)+IOIs(i+1));
    end
    as.ratio_notes{nr,1} = ratio;
    as.ratio_noteMean(nr,1) = nanmean(ratio);
    as.ratio_noteDiffs{nr,1} = diff(ratio);
    as.ratio_noteMeanDiff(nr,1) = nanmean(abs(as.ratio_noteDiffs{nr,1})); % fixed this line 2022.03.20. had been ratio, should have been ratio diff
    
    %-------------------------------------------------------------%
    % Pulse clarity
    [r, ac] = mirpulseclarity(thisWavName);%, 'Frame', 5, 's', 10, '%');
    rd = mirgetdata(r);
    %rac = mirgetdata(ac); % autocorr if we want it
    as.pulseClarity(nr,1) = rd;
    
    %-------------------------------------------------------------%
    % increment row counter
    nr = nr+1;
    
    if toplot
        % Save our figure for this stimulus
        outfname = strcat(figpath, strcat('/', date, '_', audiofiles(ifile).name(1:end-4), '.png'));
        print('-dpng', '-r300', outfname)
    end
end

%% Save raw features
outfname = fullfile(cwd, strcat('/tables/', date, '_acousticFeatures_raw.csv'));
writetable(as, outfname)

%% Make illustrative plot 
istim = 24;
figure()
plot(as.audio{istim}, 'k')

%% Create table of means for Madita
mas = table;
mas.stim = as.stim;
% Clean up mas labels
for i = 1:length(mas.stim)
    %mas.stim{i} = mas.stim{1}(2:end); % remove the d
    mas.stim_num{i} = str2double(regexp(mas.stim{i},'\d*','Match'));
    if contains(mas.stim{i}, 'M')
        mas.stim_cat{i} = 'M';
        mas.stim_cat_num(i) = 1;
    else
        mas.stim_cat{i} = 'S';
        mas.stim_cat_num(i) = 2;
    end
    
end
% now add means of all our variables
mas.intensity_mean = as.env_noteMean;
mas.intensity_meanDiff_betweenNotes = as.env_noteMeanDiff;
mas.pitch_mean = as.pitch_noteMean;
mas.pitch_meanDiff_betweenNotes = as.pitch_noteMeanDiff;
mas.timbre_mean = as.ent_noteMean;
mas.timbre_meanDiff_betweenNotes = as.ent_noteMeanDiff;
mas.IOI_mean = as.IOI_noteMean;
mas.IOI_meanDiff_betweenNotes = as.IOI_noteMeanDiff;
mas.ratio_mean = as.ratio_noteMean;
mas.ratio_meanDiff_betweenNotes = as.ratio_noteMeanDiff;
mas.pulseClarity = as.pulseClarity;
mas.spectralFlux = as.fluxMean;
mas.recurrence = as.recurr;
mas.meanLenRec = as.meanL;
mas.maxLenRec = as.maxL;
mas.determinism = as.determ;


% sort table
mas = sortrows(mas, {'stim_cat_num';'stim_num'});

% add previously calculated AMS to this table
load(strcat(cwd,'/tables/ams_peaks')); % loads into variable 'peakFreqs' 1:15 M , then 1:15 s
mas.ampModSpectrum_peak = peakFreqs;
% Commented out so Klaus can run. 

% remove unecessary column
mas.stim_cat_num = [];

% save our table
outfname = fullfile(cwd, strcat('/tables/', date, '_acousticFeatures.csv'));
writetable(mas, outfname)
