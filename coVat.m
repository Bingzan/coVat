function [mRearrangedOrig, vPermR, vPermC, vPermU, unionRC] = coVat(mOrig, bVisualise,figureName)
%coVat Implementation of co-Vat algorithm
%
% @author: Bingzan Liang
% Last update: 25/03/2014
%
    row = size(mOrig,1);
    column = size(mOrig,2);
    mDrTemp = zeros(row,row);
    mDcTemp = zeros(column,column);
    mRearrangedOrig = zeros(row, column);
    
    for i=1:row
        for j=1:row
            mDrTemp(i,j) = norm(mOrig(i,:)-mOrig(j,:));
        end
    end
    
    for i=1:column
        for j=1:column
            mDcTemp(i,j) = norm(mOrig(:,i)-mOrig(:,j));
        end
    end
    
    meanOrig = mean(mean(mOrig));
    meanDr = mean(mean(mDrTemp));
    meanDc = mean(mean(mDcTemp));
    
    lambdar = meanOrig / meanDr;
    lambdac = meanOrig / meanDc;
    
    mDr = mDrTemp% * lambdar;
    mDc = mDcTemp% * lambdac;
    
    unionRC = [[mDr mOrig];[mOrig' mDc]];
    
    [mRearrangedDr, vPermR, ~] = Vat2(mDr);
    [mRearrangedDc, vPermC, ~] = Vat2(mDc);
    [mRearrangedDu, vPermU, ~] = Vat2(unionRC);
    
    for i=1:row
        for j=1:column
            mRearrangedOrig(i,j) = mOrig(vPermR(i),vPermC(j));
        end
    end
    
    if bVisualise
       vatFigure = figure;
       colormap(gray);
       imagesc(mRearrangedOrig);
       saveas(vatFigure,figureName,'fig');
    end

end