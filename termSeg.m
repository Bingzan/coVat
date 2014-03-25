function [mTermSeg,terms] = termSeg(terms1,terms2,terms3,terms4)
%
% Create the matrix for term-segment
%
% @author: Bingzan Liang
% Last update: 25/03/2014
%
    terms = terms1;
    n = size(terms,1);
    for i = 1:size(terms2,1)
        if isempty(strmatch(terms2{i,1},terms,'exact'))
            terms{(n+1),1} = terms2{i,1};
            n = n+1;
        end
    end
    
    for i = 1:size(terms3,1)
        if isempty(strmatch(terms3{i,1},terms,'exact'))
            terms{(n+1),1} = terms3{i,1};
            n = n+1;
        end
    end
    
    for i = 1:size(terms4,1)
        if isempty(strmatch(terms4{i,1},terms,'exact'))
            terms{(n+1),1} = terms4{i,1};
            n = n+1;
        end
    end
    
    mTermSeg = zeros(size(terms,1),4);
    
    for i = 1:size(terms,1)
        if ~isempty(strmatch(terms{i,1},terms1,'exact'))
            mTermSeg(i,1) = 1;
        end
        
        if ~isempty(strmatch(terms{i,1},terms2,'exact'))
            mTermSeg(i,2) = 1;
        end
        
        if ~isempty(strmatch(terms{i,1},terms3,'exact'))
            mTermSeg(i,3) = 1;
        end
        
        if ~isempty(strmatch(terms{i,1},terms4,'exact'))
            mTermSeg(i,4) = 1;
        end
    end
    mTermSeg = 1 - mTermSeg;

end