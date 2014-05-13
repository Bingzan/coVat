function [ mSimilarity ] = similarity_cos( mAdjMat )
%SIMILARITY_COS Calculate the cosine distances among the points in the
%matrix
%
% @author: Bingzan Liang
% Last update: 13/05/2014


    mSimilarity = zeros(size(mAdjMat,1), size(mAdjMat,2));
    for i = 1:size(mAdjMat,2)
        sumRow = 0;
        for j = 1:size(mAdjMat,1)
            sumRow = sumRow + mAdjMat(j,i);
        end
        len = sqrt(sumRow);
        if len == 0
            len = size(mAdjMat,2);
        end
        for j = 1:size(mAdjMat,1)
            mAdjMat(j,i) = mAdjMat(j,i) / len;
        end
    end
    
    for i = 1:size(mAdjMat,1)
        for j = 1:size(mAdjMat,2)
            mSimilarity(i,j) = sum(mAdjMat(:,i) .* mAdjMat(:,j));
        end
    end
end

