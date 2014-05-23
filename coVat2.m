function [mRearrangedOrig1, mRearrangedOrig2, vPermR, vPermR1, vPermR2,...
    vPermC1, vPermC2, vPermU1, vPermU2, unionRC1, unionRC2] = coVat2(mOrig1, mOrig2, alpha, figureName, bVisualise)
%coVat2 Visualise each dimension according to the factor alpha
%
% @author: Bingzan Liang
% Last update: 13/05/2014
%
    %[mOrig1, mOrig2] = getData(); %termDoc, termSeg
    
    %matlabpool local 4

    %similarity
    row1 = size(mOrig1,1);
    column1 = size(mOrig1,2);
    mDr1 = zeros(row1,row1);
    mDc1 = zeros(column1,column1);
    mRearrangedOrig1 = zeros(row1, column1);
    
    %dissimilarity
    row2 = size(mOrig2,1);
    column2 = size(mOrig2,2);
    mDr2 = zeros(row2,row2);
    mDc2 = zeros(column2,column2);
    mRearrangedOrig2 = zeros(row2, column2);
    
    vPerm = cell(1,7);
    mMatrices = cell(1,7);
    
    tic
    parfor i=1:row1
        for j=1:row1
            mDr1(i,j) = norm(mOrig1(i,:)-mOrig1(j,:));
            %mDr1(i,j) = sum(mOrig1(i,:) .* mOrig1(j,:));
        end
    end
    
    parfor i=1:row2
        for j=1:row2
            mDr2(i,j) = norm(mOrig2(i,:)-mOrig2(j,:));
            %mDr2(i,j) = sum(mOrig2(i,:) .* mOrig2(j,:));
        end
    end
    
    parfor i=1:column1
        for j=1:column1
            mDc1(i,j) = norm(mOrig1(:,i)-mOrig1(:,j));
            %mDc1(i,j) = sum(mOrig1(:,i) .* mOrig1(:,j));
        end
    end
    
    parfor i=1:column2
        for j=1:column2
            mDc2(i,j) = norm(mOrig2(:,i)-mOrig2(:,j));
            %mDc2(i,j) = sum(mOrig2(:,i) .* mOrig2(:,j));
        end
    end
    
    mDr = (1-alpha)*mDr1 + alpha * mDr2;
    %Vat(mDr1,1);
    %Vat(mDr2,1);
    %Vat(mDr,1);
    
    unionRC1 = [[mDr mOrig1];[mOrig1' mDc1]];
    unionRC2 = [[mDr mOrig2];[mOrig2' mDc2]];
    
    mMatrices{1,1} = mDr;
    mMatrices{1,2} = mDr1;
    mMatrices{1,3} = mDc1;
    mMatrices{1,4} = unionRC1;
    mMatrices{1,5} = mDr2;
    mMatrices{1,6} = mDc2;
    mMatrices{1,7} = unionRC2;
    
    parfor i = 1:size(mMatrices,2)
        [~,vPerm{1,i}] = Vat2(mMatrices{1,i});
    end
    
    vPermR = vPerm{1,1};
    vPermC1 = vPerm{1,3};
    vPermC2 = vPerm{1,6};
    vPermR1 = vPerm{1,2};
    vPermU1 = vPerm{1,4};
    vPermR2 = vPerm{1,5};
    vPermU2 = vPerm{1,7};
    
    parfor i=1:row1
        for j=1:column1
            mRearrangedOrig1(i,j) = mOrig1(vPermR(i),vPermC1(j));
        end
    end
    %mRearrangedOrig1 = densify( mRearrangedOrig1,51,0.01 );
    mRearrangedOrig1 = -mRearrangedOrig1;
    mRearrangedOrig1 = densify(mRearrangedOrig1,7,151,-0.5 );
    
    
    parfor i=1:row2
        for j=1:column2
            mRearrangedOrig2(i,j) = mOrig2(vPermR(i),vPermC2(j));
        end
    end
    
    mRearrangedOrig2 = densify( mRearrangedOrig2,3,1,0.5 );
    toc
    
    if bVisualise
       vatFigure = figure;
       %colormap(hot);
       subplot 121, imagesc(mRearrangedOrig1), colormap(hot);
       freezeColors;
       subplot 122, imagesc(mRearrangedOrig2), colormap(hot);
       saveas(vatFigure,figureName,'fig');
    end

end