function [ mNewMatrix ] = densify( mRearranged,sizeAddRow,sizeAddColumn,incremental )
%DENSIFY Dense the sparse matrix
%
%Input:
%sizeAdd: the size of the density matrix (better to be an odd number)
%incremental: the incremental factor
%
%Output:
%mNewMatrix: the new matirx being densified
%
% @author: Bingzan Liang
% Last update: 2/05/2014

    mNewMatrix = mRearranged;
    row = size(mNewMatrix,1);
    column = size(mNewMatrix,2);
    halfRow = (sizeAddRow-1)/2;
    halfColumn = (sizeAddColumn-1)/2;
    [x,y,~] = find(mNewMatrix);
    for i=1:length(x)
        x1 = x(i) - halfRow;
        x2 = x(i) + halfRow;
        y1 = y(i) - halfColumn;
        y2 = y(i) + halfColumn;
                
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