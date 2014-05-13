function [ mSimilarity ] = similarity_norm( mAdjMat )
%SIMILARITY_NORM Calculate the similarity matrix usint norm-2
%   
% @author: Bingzan Liang
% Last update: 13/05/2014
%

    mSimilarity = zeros(size(mAdjMat,1), size(mAdjMat,2));
    for i = 1:size(mAdjMat,1)
        for j = 1:size(mAdjMat,2)
            mSimilarity(i,j) = norm(mAdjMat(:,i) - mAdjMat(:,j));
        end
    end
    
end

