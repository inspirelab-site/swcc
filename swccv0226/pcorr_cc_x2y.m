
function [cc_bic, corr_bic, Nh_bic, cc_aic, corr_aic, Nh_aic, corr_std, corr0, H]=pcorr_cc_x2y(V1, V2, Constraint)
    Ntimepoints=length(V2); Nh=length(V1)-Ntimepoints+1; 
    ML=zeros(Nh,1); AIC_coef=zeros(Nh,1);

    %-- Generate L_virtual Matrix from V1--BEGIN--
    L_virtual=zeros(Ntimepoints,Nh);
    for k=1:1:Ntimepoints  
    %     L_virtual(k,1:k)=V1(k:(-1):1, 1);    
        L_virtual(k,1:Nh)=V1(k+Nh-1:-1:k);
    end
    %-- Generate L_virtual Matrix from V1--END---

    %% compute all h or hps 
    H=zeros(Nh,Nh); options=optimset('Display', 'off');
    for iNh=1:Nh
        L=L_virtual(:,1:iNh);

        [h, flg] = lsqr(L,V2,1e-6,100); %unconstraint LS h    
        switch(Constraint)
            case 'all'
                H(iNh, 1:iNh)=double(h_tmp)';
                divV1V2=V2-L*h; %TO BE WRITTEN
            case 'p'  
                hp_min=min(h(h>0)); if isempty(hp_min)==1, hp_min=1e-07; end
                h_st=h; h_st(h<=0)=min(hp_min);
                hp=lsqlin(L, V2, [],[],[],[], zeros(length(h),1),[], h_st, options); %positively constrained LS h
                H(iNh, 1:iNh)=double(hp)';
                divV1V2=V2-L*hp; %TO BE WRITTEN
        end

        ML(iNh)=-Ntimepoints/2*log(2*pi*(divV1V2'*divV1V2)/(Ntimepoints-iNh))-(Ntimepoints-iNh)/2; % the maximum likelihood
        AIC_coef(iNh)=2*iNh-2*ML(iNh);

    end


    %% search for the best iNh
    [AIC_matlab, BIC_matlab]=aicbic(ML, [1:Nh], Ntimepoints);
    AICc_matlab=AIC_coef+(2*[1:Nh].*([1:Nh]+1)./(Ntimepoints-[1:Nh]-1))';

    [~, Nh_aic_matlab]=min(AIC_matlab);
    [~, Nh_aicc_matlab]=min(AICc_matlab);
    [~, Nh_bic]=min(BIC_matlab);

    if (Ntimepoints/Nh_aic_matlab)<=40
        Nh_aic=Nh_aicc_matlab;
    else
        Nh_aic=Nh_aic_matlab;
    end                

    L0=L_virtual(:,1);
    L_aic_matlab=L_virtual(:, 1:Nh_aic);
    L_bic_matlab=L_virtual(:, 1:Nh_bic);

    h_aic_matlab=squeeze(H(Nh_aic, 1:Nh_aic));h_aic_matlab=h_aic_matlab(:); %IMPULSE RESPONSE  
    V2_preAIC_matlab=L_aic_matlab*h_aic_matlab; 
    h_bic_matlab=squeeze(H(Nh_bic, 1:Nh_bic)); h_bic_matlab=h_bic_matlab(:);%IMPULSE RESPONSE  
    V2_preBIC_matlab=L_bic_matlab*h_bic_matlab; 
    h0=H(1,1); %IMPULSE RESPONSE  
    V2_pre0=L0*h0;    

    [coef0,~]=corrcoef(V2, V2_pre0);
    corr0=coef0(1,2);   

    [coef_aic_matlab,~]=corrcoef(V2, V2_preAIC_matlab);
    corr_aic=coef_aic_matlab(1,2);  

    [coef_bic_matlab,~]=corrcoef(V2, V2_preBIC_matlab);
    corr_bic=coef_bic_matlab(1,2);  

    [coef_std,~]=corrcoef(V2, V1(Nh:end));
    corr_std=coef_std(1,2);


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
        cc_aic = cap; % bit per sample 
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
        cc_bic = cap; % bit per sample 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
