function [cmap] = paired(maplength);
% from Colorbrewer
% Ilker Fer
% June 2016
cmap=[166,206,227
31,120,180
178,223,138
51,160,44
251,154,153
227,26,28
253,191,111
255,127,0]./255;

if nargin==1
    
len=length(cmap);
new=ones(maplength,3);
for i=1:3
    new(:,i)=interp1([1:len]./len,cmap(:,i),[1:maplength]./(maplength));
end
cmap=denan(new);
end  
 
