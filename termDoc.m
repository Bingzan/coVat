function [ mTermDoc ] = termDoc( freqMatrix )
%termDoc Create the term-document matrix according to the frequent matrix
%
% @author: Bingzan Liang
% Last update: 28/03/2014

    mTermDoc = zeros(max(freqMatrix(:,2)),max(freqMatrix(:,1))); %rows - num of docs; columns - num of terms

    for i=1:size(freqMatrix,1)
        mTermDoc(freqMatrix(i,2),freqMatrix(i,1)) = freqMatrix(i,3);
    end
    
    for i=1:size(mTermDoc,2)
        mTermDoc(:,i) = log(mTermDoc(:,i)+1) * log(max(freqMatrix(:,2)) / size(find(freqMatrix(:,1) == i),1)); %calculate the TF_IDF score
    end 
    
    filter = findMostFreq( freqMatrix,3 ); 
    
    mTermDoc = mTermDoc(:,(filter));
    mTermDoc = -mTermDoc;

   
end

