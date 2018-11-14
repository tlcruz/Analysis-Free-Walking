function [xmax] = GetPks(data,threshold)
if nargin < 1; error('Need data'); end;
C=size(data,2);
pp1=[data(1,:);data(1:end-1,:)];
pp2=[data(2:end,:);data(end,:)];
xmax(1:C)=struct('loc',[]);
for ch=1:C,
   if nargin ==1
     xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0)];
   else
     xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0 & data(:,ch)>threshold)];
   end
end
end

