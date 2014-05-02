function [ mNewMatrix ] = densify( mRearranged,sizeAdd,incremental )
%densify Density the sparse matrix
%
%input: 
%sizeAdd: the size of the densify matrix
%incremental: the incremental factor
%
%output;
%mNewMatrix: the new matrix being densified
%
% @author: Bingzan Liang
% Last updated: 2/5/2014

    mNewMatrix = mRearranged;
    row = size(mNewMatrix,1);
    column = size(mNewMatrix,2);
    half = (sizeAdd-1)/2;
    [x,y,~] = find(mNewMatrix);
    for i=1:length(x)
        x1 = x(i) - half;
        x2 = x(i) + half;
        y1 = y(i) - half;
        y2 = y(i) + half;
                
        if (x1 <= 0 && x2 <= row && y1 <= 0 && y2 <= column)
            x1 = 1;
            y1 = 1;
        end
        if (x1 <= 0 && x2 <= row && y1 >0 && y2 <= column)
            x1 = 1;
        end
        if (x1 >0 && x2 <= row && y1 <= 0 && y2 <= column)
            y1 = 1;
        end
        if (x1 <= 0 && x2 <= row && y1 > 0 && y2 > column)
            x1 = 1;
            y2 = column;
        end
        if (x1 > 0 && x2 <= row && y1 > 0 && y2 > column)
            y2 = column;
        end
        if (x1 > 0 && x2 > row && y1 > 0 && y2 <= column)
            x2 = row;
        end
        if (x1 > 0 && x2 > row && y1 > 0 && y2 > column)
            x2 = row;
            y2 = column;
        end
        if (x1 > 0 && x2 > row && y1 <= 0 && y2 <= column)
            x2 = row;
            y1 = 1;
        end
        mAdd = addi(x1,x2, y1, y2, incremental);
        mNewMatrix((x1:x2),(y1:y2)) = mNewMatrix((x1:x2),(y1:y2))+mAdd;
    end
end

function [mAdd] = addi(x1,x2,y1,y2,incremental)
    row = x2 - x1 + 1;
    column = y2 - y1 +1;
    mAdd = zeros(row, column);
    mAdd = mAdd + incremental;
end