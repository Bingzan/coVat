function [ mDocCat] = docCat( docs )
%docCat Create the document-category matrix
%   
% @author: Bingzan Liang
% Last update: 25/03/2014
%

    mDocCat = zeros(size(docs,1),5);
    mDocCat((1:351),1) = 1;
    mDocCat((352:606),2) = 1;
    mDocCat((607:954),3) = 1;
    mDocCat((955:1235),4) = 1;
    mDocCat((1236:1543),5) = 1;
    
    mDocCat = 1 - mDocCat;


end

