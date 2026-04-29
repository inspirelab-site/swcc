%%%%%%%%%%%%%%%%%%%pCorr_CC%%%%%%%%%%%%%%%%%%%%%%%
% Written by Jianan Jian Feb 24, 2026
% The code is adpated for channel capacity estimation from the code 
% published in Nan Xu, R Nathan Spreng, and Peter C Doerschuk. "Initial 
% validation for the estimation of resting-state fmri effective 
% connectivity by a generalization of the correlation approach." Frontiers 
% in Neuroscience, 11:271, 2017.

%% Input:
% IMG (Ntimepoints X Nnodes): Data input

% BOLDMaxlength: the maximum duration (in seconds) required for the caual effectiveness; this number has to be smaller than Ntimepoints*TR.

% TR: the time sampling rate (in seconds) for the scan 

% Constraint:   either 'p' or 'all' ;
%               'p': only positive impulse responses are considered; 
%               'all': both positive and negative impulse responses are considered

%% Output:
% CC_BIC (Nnodes X Nnodes): channel capacity result using BIC model

% Correlation_BIC (Nnodes X Nnodes): p-correlation result using BIC model

% Nh_BIC_matlab (Nnodes X Nnodes): Number of time lags estimated by BIC

% CC_AIC (Nnodes X Nnodes): channel capacity result using AIC model

% Correlation_AIC (Nnodes X Nnodes): p-correlation result using AIC model

% Nh_AIC_matlab (Nnodes X Nnodes): Number of time lags estimated by AIC

% Correlation_STD (Nnodes X Nnodes): the standard correlation result

% Correlation0 (Nnodes X Nnodes): p-Correlation result with no lags
% allowed. In any case, Correlation0 should give the same as Correlation_STD.

% H (Nnodes X Nnodes X N_hMAX X N_hMAX): Matrix which stores the coefficients of impulse reseponse, where N_hMAX=floor(BOLDMaxlength/TR). 
% squeeze(H(i,j,Nh_BIC_matlab(i,j),:)): the vector of impulse response coefficients for the information flow i-->j, which are estimated based on BIC
% squeeze(H(i,j,Nh_AIC_matlab(i,j),:)): the vector of impulse response coefficients for the information flow i-->j, which are estimated based on AIC

%Except for Correlation_STD and Correlation0, all output matrices are asymmetric

%% Example: Compute the p-correlation between timeseries x and timeseries y with postive constraint on the impulse responses.
% x=[0.4218; 0.9157; 0.7922; 0.9595; 0.6557]; 
% y=[0.9058; 0.1270; 0.9134; 0.6324; 0.0975]; 
% IMG=[x, y];
% TR=1.5; 
% BOLDMaxlength=5; % Let the longgest time required for the causal relation between x and y is 5 seconds. So, N_hMax=floor(BOLDMaxlength/TR)=3.
% Constraint='p';

% [CC_BIC, Correlation_BIC, Nh_BIC_matlab, CC_AIC, Correlation_AIC, Nh_AIC_matlab, Correlation_STD, Correlation0, H]=pCorr_CC(IMG, BOLDMaxlength, TR, Constraint);

% Correlation_STD(1,2)=-0.3185; % It computes corr(x,y)

%%information flow from x to y (x-->y)
% CC_BIC(1,2)=0.8522; %the connectivity strength determined by the channel capacity (based on BIC)
% Correlation_BIC(1,2)=-0.3185; %the connectivity strength determined by the p-Corr value (based on BIC)
% Nh_BIC_matlab(1,2)=1 %the duration of the information flow is 1*TR=1.5 secs (based on BIC)
% CC_AIC(1,2)=0.8522; %the connectivity strength determined by the channel capacity (based on AIC)
% Correlation_AIC(1,2)=-0.3185; %the connectivity strength determined by the p-Corr value (based on AIC)
% Nh_AIC_matlab(1,2)=1 %the duration of the information flow is 1*TR=1.5 secs (based on AIC)
%So, the vector of impulse response coefficients estimated both by AIC and by BIC is squeeze(H_store(1, 2, 1, :))

%%information flow from y to x (y-->x)
% CC_BIC(2,1)=2.9593; %the connectivity strength determined by the channel capacity (based on BIC)
% Correlation_BIC(2,1)=0.9007; %the connectivity strength determined by the p-Corr value (based on BIC)
% Nh_BIC_matlab(2,1)=3 %the duration of the information flow is 3*TR=4.5 secs (based on BIC)
%So, the vector of impulse response coefficients estimated by BIC is squeeze(H_store(2, 1, 3, :))
% CC_AIC(2,1)=2.7947; %the connectivity strength determined by the channel capacity (based on AIC)
% Correlation_AIC(2,1)=0.8362; %the connectivity strength determined by the p-Corr value (based on AIC)
% Nh_AIC_matlab(2,1)=2 %the duration of the information flow is 2*TR=3 secs (based on AIC)
%So, the vector of impulse response coefficients estimated by AIC is squeeze(H_store(2, 1, 2, :))


%% Program
function [CC_BIC, Correlation_BIC, Nh_BIC_matlab, CC_AIC, Correlation_AIC, Nh_AIC_matlab, Correlation_STD, Correlation0, H]=pCorr_CC(IMG, BOLDMaxlength, TR, Constraint)

N_hMAX=floor(BOLDMaxlength/TR);
N_h=1:N_hMAX;
[Ntimepoints,Nnodes]=size(IMG);

options=optimset('Display', 'off');

AIC_coef=zeros(Nnodes, Nnodes, length(N_h));
AIC_matlab=zeros(Nnodes, Nnodes, length(N_h));
AICc_matlab=zeros(Nnodes, Nnodes, length(N_h));
BIC_matlab=zeros(Nnodes, Nnodes, length(N_h));

Correlation_STD=zeros(Nnodes, Nnodes); %TO BE SAVED
Correlation0=zeros(Nnodes, Nnodes); %TO BE SAVED
Correlation_AIC=zeros(Nnodes, Nnodes); %TO BE SAVED
Correlation_BIC=zeros(Nnodes, Nnodes); %TO BE SAVED
Nh_AIC_matlab=zeros(Nnodes, Nnodes);
Nh_BIC_matlab=zeros(Nnodes, Nnodes);
CC_AIC=zeros(Nnodes, Nnodes);
CC_BIC=zeros(Nnodes, Nnodes);

ML=zeros(Nnodes, Nnodes, length(N_h));

switch(Constraint)
    case 'all'
        H_store=zeros(Nnodes, Nnodes, length(N_h), length(N_h));
    case 'p'  
        Hp_store=zeros(Nnodes, Nnodes, length(N_h), length(N_h));
end


Nh_counter=0;
for Nh=N_h
    Nh_counter=Nh_counter+1;

    for in=1:Nnodes
        V1=IMG(:,in);

        %-- Generate L_virtual Matrix from V1--BEGIN--
        L_virtual=zeros(Ntimepoints,Ntimepoints);
        for k=1:1:Ntimepoints  
            L_virtual(k,1:k)=V1(k:(-1):1, 1);
        end
        %-- Generate L_virtual Matrix from V1--END---

        L=L_virtual(:,1:Nh);
        Lpinv=pinv(L);               
        for out=1:Nnodes
            if in==out
                continue
            else
                V2=IMG(:,out);                  

                h_tmp=zeros(length(N_h),1);
                hp_tmp=zeros(length(N_h),1);
                h=Lpinv*V2; %TO BE WRITTEN
                h_tmp(1:Nh)=h;
                
                h_temp=h(1:Nh); h_st=h(1:Nh);
                h_temp(h_temp<=0)=[];
                if isempty(h_temp)==1
                    h_temp=1e-07;
                end
                h_st(h>0)=h(h>0);
                h_st(h<=0)=min(h_temp);
                hp=lsqlin(L, V2, [],[],[],[], zeros(length(h),1),[], h_st, options);  
                hp_tmp(1:Nh)=hp;                
                
                switch(Constraint)
                    case 'all'
                        H_store(in, out, Nh, :)=double(h_tmp);
                        divV1V2=V2-L*h; %TO BE WRITTEN
                    case 'p'                        
                        Hp_store(in, out, Nh, :)=double(hp_tmp);
                        divV1V2=V2-L*hp; %TO BE WRITTEN
                end
            end

            ML(in, out, Nh_counter)=-Ntimepoints/2*log(2*pi*(divV1V2'*divV1V2)/(Ntimepoints-Nh))-(Ntimepoints-Nh)/2; % the maximum likelihood
            AIC_coef(in, out, Nh_counter)=2*Nh-2*ML(in, out, Nh_counter);

        end                
        clear L Lpinv
    end

end % End of Nh variation

% save(['./H_storage/sim' sim '/subj' subj], 'H_store', 'Hp_store');
switch(Constraint)
    case 'all'
        H=H_store;
    case 'p'  
        H=Hp_store;
end

for in=1:Nnodes
    V1=IMG(:,in);

    %-- Generate L_virtual Matrix from V1--BEGIN--
    L_virtual=zeros(Ntimepoints,Ntimepoints);
    for k=1:1:Ntimepoints  
        L_virtual(k,1:k)=V1(k:(-1):1, 1);
    end
    %-- Generate L_virtual Matrix from V1--END---

    for out=1:Nnodes
        if in==out
            continue
        end
        V2=IMG(:,out);

        [AIC_matlab(in, out, :), BIC_matlab(in, out, :)]=aicbic(squeeze(ML(in,out,:)), N_h, Ntimepoints);
        AICc_matlab(in, out, :)=squeeze(AIC_coef(in, out, :))+(2*N_h.*(N_h+1)./(Ntimepoints-N_h-1))';

        [~, Nh_aic_matlab]=min(AIC_matlab(in, out, :));
        [~, Nh_aicc_matlab]=min(AICc_matlab(in, out, :));
        [~, Nh_BIC_matlab(in, out)]=min(BIC_matlab(in, out, :));

        if (Ntimepoints/Nh_aic_matlab)<=40
            Nh_AIC_matlab(in, out)=Nh_aicc_matlab;
        else
            Nh_AIC_matlab(in, out)=Nh_aic_matlab;
        end                

        L0=L_virtual(:,1);
        L_aic_matlab=L_virtual(:, 1:Nh_AIC_matlab(in, out));
        L_bic_matlab=L_virtual(:, 1:Nh_BIC_matlab(in, out));

        h_aic_matlab=squeeze(H(in,out, Nh_AIC_matlab(in, out), 1:Nh_AIC_matlab(in, out))); %IMPULSE RESPONSE  
        V2_preAIC_matlab=L_aic_matlab*h_aic_matlab; 
        h_bic_matlab=squeeze(H(in,out,Nh_BIC_matlab(in, out), 1:Nh_BIC_matlab(in, out))); %IMPULSE RESPONSE  
        V2_preBIC_matlab=L_bic_matlab*h_bic_matlab; 
        h0=H(in,out,1,1); %IMPULSE RESPONSE  
        V2_pre0=L0*h0;    

        [coef0,~]=corrcoef(V2, V2_pre0);
        Correlation0(in, out)=coef0(1,2);   

        [coef_aic_matlab,~]=corrcoef(V2, V2_preAIC_matlab);
        Correlation_AIC(in, out)=coef_aic_matlab(1,2);  

        [coef_bic_matlab,~]=corrcoef(V2, V2_preBIC_matlab);
        Correlation_BIC(in, out)=coef_bic_matlab(1,2);  

        [coef_std,~]=corrcoef(V2, V1);
        Correlation_STD(in, out)=coef_std(1,2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %AIC
        divV1V2 = V2 - V2_preAIC_matlab; % residual                                                
        v = (divV1V2' * divV1V2) / Ntimepoints; % MLE of the noise variance
        n = 512; G = fft(h_aic_matlab, n); S = v ./ abs(G) .^ 2; % EIN SPD 

        S = sort(S);
        for i = 1 : n
            lambda = n / i + mean(S(1 : i));
            if i == n || lambda <= S(i + 1)
                cap = sum(log2(lambda ./ S(1 : i))) / (2 * n);
                break;
            end
        end
        CC_AIC(in, out) = cap; % bit per sample 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %BIC
        divV1V2 = V2 - V2_preBIC_matlab; % residual                                            
        v = (divV1V2' * divV1V2) / Ntimepoints; % MLE of the variance
        n = 512; G = fft(h_bic_matlab, n); S = v ./ abs(G) .^ 2; % EIN SPD 
    
        S = sort(S);
        for i = 1 : n
            lambda = n / i + mean(S(1 : i));
            if i == n || lambda <= S(i + 1)
                cap = sum(log2(lambda ./ S(1 : i))) / (2 * n);
                break;
            end
        end
        CC_BIC(in, out) = cap; % bit per sample 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

end
