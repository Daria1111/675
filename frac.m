%clear all;
clc;
close all;
%load('gendata678new.mat')
load('/MATLAB Drive/project/3_17_0003_04_08_1.mat')

rng(42);
%data = cat(1,s1,s2);
%data = cat(1,data,s3);

chSel = [1:15];

data = data(1:15,1:200);
data = bsxfun(@minus, data, mean(data,2));
data = data(chSel,:);
chTotal = length(chSel);

order = zeros(1,chTotal);
for i = 1:chTotal
    order(i) = WT_estimator_v3(data(i,:),1);
end
baseInd = ceil(length(chSel)/2)+1;
numInp = ceil(baseInd/2);
%indTuple = {1:baseInd,setdiff(1:chTotal,1:baseInd)};
%{
for i = baseInd+1:chTotal-1
    indTuple = [indTuple;{1:i, setdiff(1:chTotal,1:i)}];
end
%}
indTuple = {[1,2,3,4,5,6,7,8,9,10],[11,12,13,14,15]};
relErr = zeros(size(indTuple,1),2);
numInp = 0;
%baseInd = 12;
%for i = 1:size(indTuple,1)
for i = 1
    relErr(i,:) = latentFracUU(data, order, numInp, indTuple{i,1}, indTuple{i,2}, 1:baseInd);
end
%}