function [ mTermDoc ] = coTermDoc(mFreqMatrix, terms)
% coTermDoc Create the matrix for term-document considering all the terms
%   mFreqmatrix: the frequent matrix
%   terms: all the terms 
%
% @author: Bingzan Liang
% last update: 25/03/2014

    mTermDoc = zeros(size(terms,1),max(mFreqMatrix(:,2)));
    
    for i=1:size(mFreqMatrix,1)
        mTermDoc(mFreqMatrix(i,1),mFreqMatrix(i,2)) = mFreqMatrix(i,3);
    end
    
    mTermDoc = -mTermDoc;

end

