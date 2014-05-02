function [mRearrangedOrig1, mRearrangedOrig2, vPermR, vPermR1, vPermR2, vPermC1, vPermC2, vPermU1, vPermU2, unionRC1, unionRC2] = coVat2(alpha, figureName, bVisualise)
%coVat2 Visualise each dimension according to the factor alpha
%
% @author: Bingzan Liang
% Last update: 1/05/2014
%
    [mOrig1, mOrig2] = getData(201,-0.01); %termDoc, termSeg

    row1 = size(mOrig1,1);
    column1 = size(mOrig1,2);
    mDrTemp1 = zeros(row1,row1);
    mDcTemp1 = zeros(column1,column1);
    mRearrangedOrig1 = zeros(row1, column1);
    
    row2 = size(mOrig2,1);
    column2 = size(mOrig2,2);
    mDrTemp2 = zeros(row2,row2);
    mDcTemp2 = zeros(column2,column2);
    mRearrangedOrig2 = zeros(row2, column2);
    
    for i=1:row1
        for j=1:row1
            mDrTemp1(i,j) = norm(mOrig1(i,:)-mOrig1(j,:));
        end
    end
    
    for i=1:row2
        for j=1:row2
            mDrTemp2(i,j) = norm(mOrig2(i,:)-mOrig2(j,:));
        end
    end
    
    for i=1:column1
        for j=1:column1
            mDcTemp1(i,j) = norm(mOrig1(:,i)-mOrig1(:,j));
        end
    end
    
    for i=1:column2
        for j=1:column2
            mDcTemp2(i,j) = norm(mOrig2(:,i)-mOrig2(:,j));
        end
    end
    
    %meanOrig = mean(mean(mOrig));
    %meanDr = mean(mean(mDrTemp));
    %meanDc = mean(mean(mDcTemp));
    
    %lambdar = meanOrig / meanDr;
    %lambdac = meanOrig / meanDc;
    
    mDr1 = mDrTemp1; % * lambdar;
    mDc1 = mDcTemp1; % * lambdac;
    mDr2 = mDrTemp2;
    mDc2 = mDcTemp2;
    
    mDr = (1-alpha)*mDr1 + alpha * mDr2;
    %Vat(mDr1,1);
    %Vat(mDr2,1);
    %Vat(mDr,1);
    
    unionRC1 = [[mDr mOrig1];[mOrig1' mDc1]];
    unionRC2 = [[mDr mOrig2];[mOrig2' mDc2]];
    
    [~, vPermR, ~] = Vat2(mDr);
    
    [~, vPermR1, ~] = Vat2(mDr1);
    [~, vPermC1, ~] = Vat2(mDc1);
    [~, vPermU1, ~] = Vat2(unionRC1);
    
    [~, vPermR2, ~] = Vat2(mDr2);
    [~, vPermC2, ~] = Vat2(mDc2);
    [~, vPermU2, ~] = Vat2(unionRC2);
    
    for i=1:row1
        for j=1:column1
            mRearrangedOrig1(i,j) = mOrig1(vPermR(i),vPermC1(j));
        end
    end
    mRearrangedOrig1 = sparse(mRearrangedOrig1);
    
    for i=1:row2
        for j=1:column2
            mRearrangedOrig2(i,j) = mOrig2(vPermR(i),vPermC2(j));
        end
    end
    %mRearrangedOrig2 = sparse(mRearrangedOrig2);
    
    if bVisualise
       vatFigure = figure;
       colormap(gray);
       subplot 121, spy(mRearrangedOrig1);
       subplot 122, imagesc(mRearrangedOrig2);
       saveas(vatFigure,figureName,'fig');
    end

end