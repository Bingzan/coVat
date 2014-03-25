function [mTermDoc, mDocCat] = getData()
% getData Read all related files and create the matrix
%
% @author: Bingzan Liang
% Last update: 25/03/2014

    %terms used in term-seg
    %terms1 = textread('bbc_seg1of4.terms','%s','delimiter','\n'); 
    %terms2 = textread('bbc_seg2of4.terms','%s','delimiter','\n'); 
    %terms3 = textread('bbc_seg3of4.terms','%s','delimiter','\n'); 
    %terms4 = textread('bbc_seg4of4.terms','%s','delimiter','\n'); 
    %[mTermSeg,terms] = termSeg(terms1,terms2,terms3,terms4);
    
    docs = textread('bbc_seg1of4.docs','%s','delimiter','\n');
    [ mDocCat] = docCat( docs );
    
    mFreqMatrix = load('bbc_seg1of4.mtx');
    mTermDoc = termDoc( mFreqMatrix );

end