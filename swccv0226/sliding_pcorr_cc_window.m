%modified version of sliding window program calculates sliding window
%p-correlation and channel capacity between two time courses.  used to 
%examine rois from a network. Sliding tc1 with memeory
% last modified by Jianan Jian on 02242026

function [cc_bic, corr_bic, Nh_bic, cc_aic, corr_aic, Nh_aic, corr_std, corr0, H, nmaps]=sliding_pcorr_cc_window(tc1,tc2, windowsize, overlap, maxLag, Constraint)
    Ntimepoints=length(tc1); tc1f=[zeros(maxLag-1,1); tc1(:)]; 
    nmaps=0; H=[]; tic
    for i=1:(windowsize-overlap):(Ntimepoints-maxLag-1)
       if (i+windowsize-1)<=Ntimepoints
            nmaps=nmaps+1;
            V1=tc1f(i:i+maxLag+windowsize-2); V2=tc2(i:i+windowsize-1);
            [cc_bic(nmaps), corr_bic(nmaps), Nh_bic(nmaps), cc_aic(nmaps), corr_aic(nmaps), Nh_aic(nmaps), corr_std(nmaps), corr0(nmaps), h]=pcorr_cc_x2y(V1, V2, Constraint);
            H=cat(3,H,h);  
       else
            nmaps=nmaps+1; windowsize2=Ntimepoints+1-i;
            V1=tc1f(i:i+maxLag+windowsize2-2); V2=tc2(i:i+windowsize2-1);
            [cc_bic(nmaps), corr_bic(nmaps), Nh_bic(nmaps), cc_aic(nmaps), corr_aic(nmaps), Nh_aic(nmaps), corr_std(nmaps), corr0(nmaps), h]=pcorr_cc_x2y(V1, V2, Constraint);
            H=cat(3,H,h);                  
       end
%        disp(['sliding window #' num2str(nmaps) ', time: ' num2str(toc)]);
    end
    
disp(['complete #' num2str(nmaps) ' sliding windows, time: ' num2str(toc)]);
end