



%studyPath      = '/path/to/study'
%Set this to the specific directory your study is. 
studyPath = '.'

%The EEGLAB "study" file.  Set this to whatever your study is.
%Default here for the EEGLAb "animal study dataset available from:
%https://sccn.ucsd.edu/eeglab/download/animal_study.zip

studyFilename = 'animal.study'

if ~exist("studyFilename",'file')
    error("Can't find study")
end


%Load an eeglab study. 
[STUDY ALLEEG] = pop_loadstudy('filename', studyFilename, 'filepath', studyPath);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

%Use the electrodes defined in the study to create  custom template
customTemplate = createCustomTemplates(EEG(1).chanlocs,EEG(1).ref)

%Can also save a mat file and reload these instead of fitting every time. 

%Read the erp data for the study
%NOTE: These need to have been precomputed using EEGLAB.  THey only get red
%in not computed. 
[STUDY, erpdata, times, setinds, cinds] = std_readdata(STUDY, ALLEEG, 'channels', { ALLEEG(1).chanlocs(:).labels },'datatype','erp');


%Average across participants. 
%Note: Need to transpose EEGlab data so that electrodes are the 
%first dimension

averageErpData = cell(size(erpdata));


for iCond = 1:length(erpdata)
    %Average across participants
    averageDataSingleCondition = mean(erpdata{iCond},3);
    
    %Set electrodes as first dimesion
    averageErpData{iCond} = permute(averageDataSingleCondition,[2 1]);   

end


%Submit average data to template fit routine:
templateFitData = fitEEGTemplates(averageErpData, customTemplate.weights);




%Plot every visual area ROI in a different sublot
%with 2 columns for the plots
%Plot every condition as a different color

absMaxFitData = max(max(abs(cell2mat(templateFitData))));
figure()

for iCond = 1:numel(templateFitData)

    numRois=size(customTemplate.listROIs,2)
    
    for iRoi=1:numRois
        subplot(ceil(numRois/2),2,iRoi)
        hold on;
        plot(times,templateFitData{iCond}(iRoi,:))
        ylim( [-absMaxFitData, absMaxFitData])
        title(customTemplate.listROIs{iRoi})
            
    end
end


