%simplified version of sliding window program calculates sliding window
%correlation between two time courses.  used to examine rois from a
%network
%needs tc1, tc2 in memory
% modified by Nan Xu

function [corrmap]=sliding_corr_window_simple(tc1,tc2, windowsize, overlap, endlength)
nmaps=1;
% windowsize= 100;
% overlap = 99;
DimTime=length(tc1);

for i=1:(windowsize-overlap):(DimTime-endlength-1)
               if ((i+windowsize-1)<=DimTime)
                    reftc=tc1(i:i+windowsize-1);
                    imgtc=tc2(i:i+windowsize-1);
                    temp = corrcoef(reftc,imgtc);
                    corrmap(1,nmaps)=temp(2,1);
                    nmaps=nmaps+1;
                else
                    reftc=tc1(i:(DimTime-endlength));
                    imgtc=tc2(i:(DimTime-endlength));
                    temp = corrcoef(reftc,imgtc);
                    corrmap(1,nmaps)=temp(2,1);
                    nmaps=nmaps+1;
               end
end