function [mRearrangedOrig, vPermR, vPermC, vPermU, unionRC] = coVat(mOrig, bVisualise,figureName)
%coVat Implementation of co-Vat algorithm
%
% @author: Bingzan Liang
% Last update: 13/05/2014
%
    %[mOrig,~,~] = svds(mOrig,1000);

    %matlabpool local 4
    
    row = size(mOrig,1);
    column = size(mOrig,2);
    mDr = zeros(row,row);
    mDc = zeros(column,column);
    mRearrangedOrig = zeros(row, column);
    
    vPerm = cell(1,3);
    mRearranged = cell(1,3);
    mMatrices = cell(1,3);
    
    tic
    parfor i=1:row
    %for i = 1:row
        for j=1:row
            mDr(i,j) = norm(mOrig(i,:)-mOrig(j,:));
            %mDr(i,j) = sum(mOrig(i,:) .* mOrig(j,:));
        end
    end
    
    parfor i=1:column
    %for i=1:column
        for j=1:column
            mDc(i,j) = norm(mOrig(:,i)-mOrig(:,j));
            %mDc(i,j) = sum(mOrig(:,i) .* mOrig(:,j));
        end
    end
    
    unionRC = [[mDr mOrig];[mOrig' mDc]];
    
    mMatrices{1,1} = mDr;
    mMatrices{1,2} = mDc;
    mMatrices{1,3} = unionRC;
    
    parfor i = 1:size(mMatrices,2)
    %for i = 1:size(mMatrices,2)
        [mRearranged{1,i}, vPerm{1,i}, ~] = Vat2(mMatrices{1,i});
    end
    
    %[mRearrangedDr, vPermR, ~] = Vat2(mDr);
    %[mRearrangedDc, vPermC, ~] = Vat2(mDc);
    %[mRearrangedDu, vPermU, ~] = Vat2(unionRC);
    %[mRearrangedDr, vPermR, ~, ~,~] = Vat(mDr, 1);
    %[mRearrangedDc, vPermC, ~, ~,~] = Vat(mDc, 1);
    
    vPermR = vPerm{1,1};
    vPermC = vPerm{1,2};
    vPermU = vPerm{1,3};
    
    parfor i=1:row
    %for i=1:row
        for j=1:column
            mRearrangedOrig(i,j) = mOrig(vPermR(i),vPermC(j));
        end
    end
    toc
    
    %mRearrangedOrig = 1 - mRearrangedOrig;
    %mRearrangedOrig = densify(mRearrangedOrig,3,3,1 );
    mRearrangedOrig = -mRearrangedOrig;
    %mRearrangedOrig = 1-mRearrangedOrig;
    
    
    if bVisualise
       vatFigure = figure;
       colormap(hot);
       imagesc(mRearrangedOrig);
       saveas(vatFigure,figureName,'fig');
    end

end