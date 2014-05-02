function [mTermDoc, mDocCat] = getData(sizeAdd, incremental)
% getData Read all related files and create the matrix
%
% @author: Bingzan Liang
% Last update: 2/05/2014
    
    docs = textread('bbc_seg1of4.docs','%s','delimiter','\n');
    [ mDocCat] = docCat( docs );
    
    mFreqMatrix = load('bbc_seg1of4.mtx');
    mTermDoc = termDoc( mFreqMatrix );
    mTermDoc = densify( mTermDoc,sizeAdd,incremental);

end