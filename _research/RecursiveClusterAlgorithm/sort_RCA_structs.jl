%% SortRCAStructs.m
%
%Matt Reeves 01/07/2014
%
%
%Puts the RCA output into a more useful form

function [PSTATS,NSTATS] = SortRCAStructs(PCLUSTERS,NCLUSTERS)

%% Sort out the positive clusters 
names = fieldnames(PCLUSTERS);
numclusters = length(names);

kappac = zeros(numclusters,1);
RC = zeros(numclusters,2);
radii = zeros(numclusters,1);
for ii = 1:numclusters
  string1 = ['cluster' num2str(ii) '.positions'];
  eval(['temp = PCLUSTERS.' string1 ';']);
  temp = size(temp);
  kappac(ii) = temp(1);
  
  string2 = ['cluster' num2str(ii) '.RC'];
  eval(['temp = PCLUSTERS.' string2 ';']);
  RC(ii,:) = temp;
    
  string3 = ['cluster' num2str(ii) '.Radius'];
  eval(['temp = PCLUSTERS.' string3 ';']);
  radii(ii) = temp;
end

[kappac,I] = sort(kappac,'descend');
RC(:,1) = RC(I,1);
RC(:,2) = RC(I,2);
radii = radii(I);

%Charge, XC, YC, radius, clusternumber
PSTATS = [kappac RC radii I];

%% Sort out the negative clusters 
names = fieldnames(NCLUSTERS);
numclusters = length(names);

kappac = zeros(numclusters,1);
RC = zeros(numclusters,2);
radii = zeros(numclusters,1);
for ii = 1:numclusters
  string1 = ['cluster' num2str(ii) '.positions'];
  eval(['temp = NCLUSTERS.' string1 ';']);
  temp = size(temp);
  kappac(ii) = temp(1);
  
  string2 = ['cluster' num2str(ii) '.RC'];
  eval(['temp = NCLUSTERS.' string2 ';']);
  RC(ii,:) = temp;
    
  string3 = ['cluster' num2str(ii) '.Radius'];
  eval(['temp = NCLUSTERS.' string3 ';' ]);
  radii(ii) = temp;
end

[kappac,I] = sort(kappac,'descend');
RC(:,1) = RC(I,1);
RC(:,2) = RC(I,2);
radii = radii(I);

%Charge, XC, YC, radius, clusternumber
NSTATS = [kappac RC radii I];

return


