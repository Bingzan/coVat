function [ filter ] = findMostFreq( freq,num )
%findMostFreq Remove the too frequent terms (less useful)
%
% @author: Bingzan Liang
% Last update: 28/3/2014
%
    
    doc_freq = zeros(max(freq(:,1)),1);
    for i=1:size(doc_freq,1)
        doc_freq(i,1) = size(find(freq(:,1) == i),1);
    end
    
    filter = find(doc_freq(:,1)<=num);

end

