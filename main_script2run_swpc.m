%script to run all sliding window correlation on all ROIS for network analysis
%must first load in time course set
clear; close all; clc
%% Key parameter setting up
Constraint='p'; % please also try:  Constraint='all';
windowsize_duration=50; %sec
TR=.5; windowsize=floor(windowsize_duration/TR/2)*2;
overlap=windowsize-1;
maxDur=15; maxLag=floor(maxDur/TR); 
%% Loading files and intialization

addpath('./swccv0226/')
fileIN='data.mat'; load(fileIN);
[ntimepoints, nsubj]=size(BOLD_ROI1_array); nroi=2; st=1; 
if ~exist('./output','dir'), mkdir('./output/'); end
SubSelect=1:nsubj;
for isbj=SubSelect
    fileOUT=['./output/subj' num2str(isbj) '.mat'];
    
    slpcorr_length=ceil((ntimepoints-maxLag-1)/(windowsize-overlap));
    bold=[squeeze(BOLD_ROI1_array(st:ntimepoints+st-1,isbj))'; squeeze(BOLD_ROI2_array(st:ntimepoints+st-1,isbj))'];   
    x=zeros(nroi, nroi, slpcorr_length); slcc_AIC=x; slcc_BIC=x; slpcorr_AIC=x; slpcorr_BIC=x; slNh_BIC=x; slNh_AIC=x; slpcorr_STD=x; slpcorr1=x; slcorr=x; 
    
    tic
    for in = 1:nroi
        disp(['--------------- (in, time):  ' num2str(in) ', ' num2str(toc) '-----------------']);
%         parfor out=in+1:nroi     
        for out=in+1:nroi     
            if in==out, continue; end
            tc1=bold(in,:)';  tc2=bold(out,:)';

%             %%% sliding window correlation
%             disp(['---------------slcorr (in, out):  ' num2str(in) ', ' num2str(out) '-----------------']);
%             slcorr(in,out,:)= sliding_corr_window_simple(tc1,tc2, windowsize, overlap, maxLag);
            
            %%% sliding window p-corr (modified version)
            disp(['---------------slpcorr1 (in, out):  ' num2str(in) ', ' num2str(out) '-----------------']);
            [slcc_BIC(in,out,:), slpcorr_BIC(in,out,:), slNh_BIC(in,out,:), slcc_AIC(in,out,:), slpcorr_AIC(in,out,:), slNh_AIC(in,out,:), slpcorr_STD(in,out,:), slpcorr1(in,out,:), ~, slpcorr_length1]=sliding_pcorr_cc_window(tc1,tc2, windowsize, overlap, maxLag, Constraint);
            disp(['---------------slpcorr1 (out, in):  ' num2str(in) ', ' num2str(out) '-----------------']);
            [slcc_BIC(in,out,:), slpcorr_BIC(out,in,:), slNh_BIC(out,in,:), slcc_AIC(in,out,:), slpcorr_AIC(out,in,:), slNh_AIC(out,in,:), slpcorr_STD(out,in,:), slpcorr1(out,in,:), ~, ~]=sliding_pcorr_cc_window(tc2,tc1, windowsize, overlap, maxLag, Constraint);
        end        
    end
    save(fileOUT);
    disp(['--------------- FINAL TIME:  ' num2str(toc) '-----------------']);
end

