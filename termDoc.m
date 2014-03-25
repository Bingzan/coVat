function [ mTermDoc ] = termDoc( freqMatrix )
%termDoc Create the term-document matrix according to the frequent matrix
%
% @author: Bingzan Liang
% Last update: 25/03/2014

    mTermDocTemp = zeros(max(freqMatrix(:,2)),max(freqMatrix(:,1))); %rows - num of docs; columns - num of terms

    for i=1:size(freqMatrix,1)
        mTermDocTemp(freqMatrix(i,2),freqMatrix(i,1)) = freqMatrix(i,3);
    end

    mTermDocTemp = -mTermDocTemp; % transform the similatiry matrix to dissimilarity matrix

    mTermDoc = mTermDocTemp;

%sample = datasample(1:4659,2000,'Replace',false);
%mTermDoc = zeros(2000,1543);
%for i = 1:2000
%    mTermDoc(i,:) = mTermDocTemp(sample(i),:);
%end

end

