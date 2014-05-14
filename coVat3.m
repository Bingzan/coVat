function [mRearrangedOrig1, mRearrangedOrig2, mRearrangedOrig3, ...
    vPermR, vPermR1, vPermR2, vPermR3, vPermC1, vPermC2, vPermC3,...
    vPermU1, vPermU2, vPermU3, unionRC1, unionRC2, unionRC3] = coVat3(mOrig1, mOrig2, mOrig3, alpha1,alpha2, figureName, bVisualise)
%COVAT3 Visualise each dimension according to the factor alpha
%
% @author: Bingzan Liang
% Last update: 14/05/2014
%
    %[mOrig1, mOrig2] = getData(); %termDoc, termSeg

    row1 = size(mOrig1,1);
    column1 = size(mOrig1,2);
    mDr1 = zeros(row1,row1);
    mDc1 = zeros(column1,column1);
    mRearrangedOrig1 = zeros(row1, column1);
    
    row2 = size(mOrig2,1);
    column2 = size(mOrig2,2);
    mDr2 = zeros(row2,row2);
    mDc2 = zeros(column2,column2);
    mRearrangedOrig2 = zeros(row2, column2);
    
    row3 = size(mOrig3,1);
    column3 = size(mOrig3,2);
    mDr3 = zeros(row3,row3);
    mDc3 = zeros(column3,column3);
    mRearrangedOrig3 = zeros(row3, column3);
    
    for i=1:row1
        for j=1:row1
            mDr1(i,j) = norm(mOrig1(i,:)-mOrig1(j,:));
            %mDr1(i,j) = sum(mOrig1(i,:) .* mOrig1(j,:));
        end
    end
    
    for i=1:row2
        for j=1:row2
            mDr2(i,j) = norm(mOrig2(i,:)-mOrig2(j,:));
            %mDr2(i,j) = sum(mOrig2(i,:) .* mOrig2(j,:));
        end
    end
    
    for i=1:row3
        for j=1:row3
            mDr3(i,j) = norm(mOrig3(i,:)-mOrig3(j,:));
            %mDr2(i,j) = sum(mOrig2(i,:) .* mOrig2(j,:));
        end
    end
    
    for i=1:column1
        for j=1:column1
            mDc1(i,j) = norm(mOrig1(:,i)-mOrig1(:,j));
            %mDc1(i,j) = sum(mOrig1(:,i) .* mOrig1(:,j));
        end
    end
    
    for i=1:column2
        for j=1:column2
            mDc2(i,j) = norm(mOrig2(:,i)-mOrig2(:,j));
            %mDc2(i,j) = sum(mOrig2(:,i) .* mOrig2(:,j));
        end
    end
    
    for i=1:column3
        for j=1:column3
            mDc3(i,j) = norm(mOrig3(:,i)-mOrig3(:,j));
            %mDc2(i,j) = sum(mOrig2(:,i) .* mOrig2(:,j));
        end
    end
    
    mDr = alpha1*mDr1 + alpha2 * mDr2 + (1-alpha1-alpha2) * mDr3;
    
    unionRC1 = [[mDr mOrig1];[mOrig1' mDc1]];
    unionRC2 = [[mDr mOrig2];[mOrig2' mDc2]];
    unionRC3 = [[mDr mOrig3];[mOrig3' mDc3]];
    
    [~, vPermR, ~] = Vat2(mDr);
    
    [~, vPermR1, ~] = Vat2(mDr1);
    [~, vPermC1, ~] = Vat2(mDc1);
    [~, vPermU1, ~] = Vat2(unionRC1);
    
    [~, vPermR2, ~] = Vat2(mDr2);
    [~, vPermC2, ~] = Vat2(mDc2);
    [~, vPermU2, ~] = Vat2(unionRC2);
    
    [~, vPermR3, ~] = Vat2(mDr3);
    [~, vPermC3, ~] = Vat2(mDc3);
    [~, vPermU3, ~] = Vat2(unionRC3);
    
    for i=1:row1
        for j=1:column1
            mRearrangedOrig1(i,j) = mOrig1(vPermR(i),vPermC1(j));
        end
    end
    %mRearrangedOrig1 = densify( mRearrangedOrig1,51,0.01 );
    
    
    for i=1:row2
        for j=1:column2
            mRearrangedOrig2(i,j) = mOrig2(vPermR(i),vPermC2(j));
        end
    end
    %mRearrangedOrig2 = densify( mRearrangedOrig2,51,0.01 );
    
    for i=1:row3
        for j=1:column3
            mRearrangedOrig3(i,j) = mOrig3(vPermR(i),vPermC3(j));
        end
    end
    
    if bVisualise
       vatFigure = figure;
       colormap(hot);
       subplot 131, imagesc(mRearrangedOrig1);
       subplot 132, imagesc(mRearrangedOrig2);
       subplot 133, imagesc(mRearrangedOrig3);
       saveas(vatFigure,figureName,'fig');
    end

end