function [mRearrangedOrig, vPermR, vPermC, vPermU, unionRC] = coVat(mOrig, bVisualise,figureName)
%coVat Implementation of co-Vat algorithm
%
% @author: Bingzan Liang
% Last update: 13/05/2014
%
    %[mOrig,~,~] = svds(mOrig,1000);

    row = size(mOrig,1);
    column = size(mOrig,2);
    mDr = zeros(row,row);
    mDc = zeros(column,column);
    mRearrangedOrig = zeros(row, column);
    
    for i=1:row
        for j=1:row
            mDr(i,j) = norm(mOrig(i,:)-mOrig(j,:));
            %mDr(i,j) = sum(mOrig(i,:) .* mOrig(j,:));
        end
    end
    
    for i=1:column
        for j=1:column
            mDc(i,j) = norm(mOrig(:,i)-mOrig(:,j));
            %mDc(i,j) = sum(mOrig(:,i) .* mOrig(:,j));
        end
    end
    
    unionRC = [[mDr mOrig];[mOrig' mDc]];
    
    [mRearrangedDr, vPermR, ~] = Vat2(mDr);
    [mRearrangedDc, vPermC, ~] = Vat2(mDc);
    [mRearrangedDu, vPermU, ~] = Vat2(unionRC);
    
    for i=1:row
        for j=1:column
            mRearrangedOrig(i,j) = mOrig(vPermR(i),vPermC(j));
        end
    end
    
    %mRearrangedOrig = densify( mRearrangedOrig,51,0.001 );
    
    if bVisualise
       vatFigure = figure;
       colormap(gray);
       imagesc(mRearrangedOrig);
       saveas(vatFigure,figureName,'fig');
    end

end