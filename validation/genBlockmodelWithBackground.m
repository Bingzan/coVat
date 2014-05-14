function [cmAdjMat, vVertPos, cmImageGraph, cmShuffledAdjMat, vShuffledVertPos] =...
        genBlockmodelWithBackground_rec(sPosDist, posNum, numRow, numColumn, unitGraphSize, bShuffle, vBackgroundProp, sBlockmodelType,...
            varargin)
%
% Generate blockmodels, using a stocastic blockmodel (generative) model.
% Compared to genBlockmodel, this version does not need specification of sparse
% blocks (it assumes 0 density initially).  Instead, it follows Newman's
% approach and introduces a background, erdos renyi model, whose contribution is
% parameterised.
%
% It also allows generation of different predefined structure (like the newman
% version).  At the moment, these are: community, core-periphery, hierarchy,
% bipartite and random (as before).  
%
% @author Jeffrey Chan
% @Date 1/2013
%
% Modified by: Bingzan Liang
% Last update: 14/5/2014
%
% Input:
% sPosDist - probability distribution of the sizes of each position/role
% posNum - number of roles to generate
% numRow - number of row blocks (should be the same with posNumb)
% numColumn - number of column blocks
% unitGraphSize - unit number of points in each block
% bShuffle - whether to shuffle the generated adjacency matrix
% vBackgroundProp - vector of background portions to generate for each generated graph
% sBlockmodelType - the type of blockmodel to generate (community, corePeriphery, hierarchy, bipartite, random).
%
% Output:
% mImageGraph - image matrix, each entry being 0-1, specifying the
%   probability (i.e., the density of each block) of an edge existing.  The
%   mImageGraph must be a square matrix and the number of rows (columns) = size(vRoleProb,2).
% mAdjMat - (unshuffled) generated adajency matrix.
% vVertRole - mapping of vertices to roles.
% mShuffledAdjMat - shuffled adajency matrix.
% vShuffledVertRole - shuffled vertices to roles mapping.
%

    mystream = RandStream('mt19937ar','Seed',sum(100*clock));
    RandStream.setGlobalStream(mystream);

    %graphSize = max(numRow*unitGraphSize, numColumn*unitGraphSize);
    optargin = size(varargin,2);
    stdargin = nargin - optargin;

    if stdargin < 8
        disp(' genBlockmodelWithBackground_rec(sPosDist, posNum, numRow, numColumn, unitGraphSize, bShuffle, vBackgroundProp, sBlockmodelType, varargin)');
        %exit;
    end

    currIndex = 1;
    % file name
    bFileOutput = false;

    % if image matrix missing, we generate it
    if optargin > 0
        % first arguement will determine what should be done
        % generate image matrix
        % one image matrix will be generated per entry of vMeanDense
        vMeanDense = varargin{1};
        vVarDense = varargin{2};
        cOptional = [];
        if strcmp(sBlockmodelType, 'random')
            probDense = varargin{3};            
            cOptional = [probDense];
            currIndex = currIndex + 1;
        elseif strcmp(sBlockmodelType, 'corePeriphery')
            stepSize = varargin{3};
            cOptional = [stepSize];
            currIndex = currIndex + 1;
        end
        currIndex = currIndex + 2;
        cmImageGraphPlanted = genImageGraphPlanted(sBlockmodelType, posNum,numRow, numColumn, vMeanDense, vVarDense, 0.0, 0.0, cOptional);
        
        % more command line arguments
        if optargin >= currIndex
            if strcmp(varargin{currIndex}, 'file')
                bFileOutput = true;
                sAdjFilename = varargin{currIndex+1};
                sPosFilename = varargin{currIndex+2};
                sImageFilename = varargin{currIndex+3};
            end
        end
    end


    % generate role sizes
    [vPosNumRow, vPosNumColumn, vVertPos] = genPositionMembership(sPosDist, posNum, numRow, numColumn, unitGraphSize);

    assert(size(vPosNumRow,2) == size(cmImageGraphPlanted{1},1));
    assert(size(vPosNumColumn,2) == size(cmImageGraphPlanted{1},2));
    %assert(size(cmImageGraphPlanted{1},1) == size(cmImageGraphPlanted{1},2)); % square matrix




    % generate background image graph
    cmImageGraphBackground = genBackgroundImage(cmImageGraphPlanted, vPosNumRow, vPosNumColumn, numRow, numColumn, unitGraphSize);

    % compute mImageGraph
    cmImageGraph = genImageGraph(cmImageGraphPlanted, cmImageGraphBackground, vBackgroundProp);
   
    % generate the graph
    cmAdjMat = genGraph(vPosNumRow, vPosNumColumn, cmImageGraph, numRow, numColumn, unitGraphSize);
    
    
    [vShuffledVertPos, cmShuffledAdjMat] = shuffle(bShuffle, vVertPos, cmAdjMat, numRow, numColumn, unitGraphSize);
    

    % if file output, we assume it is a script and function will be terminated
    % after all processing has been finished
    if bFileOutput
        dlmwrite(sAdjFilename, mShuffledAdjMat, 'precision', '%d');
    
        % loop through the partitions and write out to file (note it appends)
        for p = 1 : posNum
            vIndices = find(vShuffledVertPos == p);
            dlmwrite(sPosFilename, vIndices - 1, '-append'); 
        end
    
    
        dlmwrite(sImageFilename, mImageGraph, 'precision', '%.2f');
       

    end

end % end of function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cmImageGraphPlanted] = genImageGraphPlanted(sBlockmodelType, posNum, numRow, numColumn, vMeanDense, vStdDense, meanSparse, stdSparse, vOptional)
%
% Generate image matrix.  Each block is either dense or sparse, depending
% on probDense (1 - probDense) respectively.  The density of each block is
% then generated by a gaussian(meanDense, stdDense) for dense (and for
% sparse, we use meanSparse and varSparse).
%
% sBlockmodelType - see description for main function.
% posNum - number of blocks
% meanDense - 
%
% vOptional (vector of extra double parameters):
% For 'random' blockmodelType
%   probDense - probability that a block is dense
% For 'coreperiphery' blockmodelType
%   stepSize - This controls the "skewness" of the core-periphery; i.e., how
%   fast does the dense block rate drop as we go further from the core
%   position?  Set to 0 for all blocks to be dense, 1 for a "typical" core
%   periphery, and higher for greater skewness.
%  
%

    
    switch sBlockmodelType
        case 'community'
            % generate dense diagonal blocks
            mDense = logical(diag(ones(1,posNum), 0));
            cmImageGraphPlanted = genBlocks(mDense, vMeanDense, vStdDense, meanSparse, stdSparse);     
        case 'corePeriphery'
            % TODO, this will generate a uniform block distributed
            % core-perhipery
            % It might be better to use newman's generation for this type of
            % structure, as it will generate graphs that are closer to reality
            stepSize = vOptional(1);
            mDense = false(posNum, posNum);
            for r = 1 : posNum
                for c = 1 : posNum - (r-1)*stepSize
                    mDense(r,c) = true;
                    mDense(c,r) = true;
                end
            end

            cmImageGraphPlanted = genBlocks(mDense, vMeanDense, vStdDense, meanSparse, stdSparse);
        case 'hierarchy'
            mDense = false(posNum, posNum);
            
            % this will generate a single 0 image matrix if posNum = 1, as
            % desired.  Otherwise generate an off-diagonal hierarchical image
            % matrix.  Note that for posNum = 2, this generates a biparatite
            % structure.
            for p = 1 : posNum-1
                mDense(p, p+1) = true;
                mDense(p+1, p) = true;
            end
            cmImageGraphPlanted = genBlocks(mDense, vMeanDense, vStdDense, meanSparse, stdSparse);
        case 'downHierarchy'
            % Hierarchy where a flow goes from top position to next one
            mDense = false(posNum, posNum);
            
            % see notes for hierarchy generation to understand the logical
            % behind this
            for p = 1 : posNum-1
                mDense(p, p+1) = true;
            end
            
            cmImageGraphPlanted = genBlocks(mDense, vMeanDense, vStdDense, meanSparse, stdSparse);
        case 'upHierarchy'
            probDense = vOptional(1);
            [cmImageGraphPlanted] = genImageGraphPlanted('downHierarchy', probDense, posNum, meanDense, stdDense, meanSparse, stdSparse);
            % upHierarchy is the tranpose of downHierarchy
            for m = 1 : size(cmImageGraphPlanted,2)
                cmImageGraphPlanted{m} = cmImageGraphPlanted{m}';
            end
        case 'random'
            probDense = vOptional(1);
            % generate dense blocks according to rand and probDense
            mDense = logical(rand(numRow, numColumn) < probDense);
            cmImageGraphPlanted = genBlocks(mDense, vMeanDense, vStdDense, meanSparse, stdSparse);
%             for r = 1 : posNum
%                 for c = 1 : posNum
%                     if mDense(r,c) < probDense
%                         density = normrnd(meanDense, stdDense);
%                         if density >= 1.0, density = 1.0; end
%                         mImageGraphPlanted(r,c) = density;
%                     else
%                         density = normrnd(meanSparse, stdSparse);
%                         if density <= 0.0, density = 0.0; end            
%                         mImageGraphPlanted(r,c) = density;
%                     end
%                 end
%             end            
        otherwise
            disp(usage());  
            err = MException('genImageGraphPlanted:InvalidBlockType', 'Invalid blockmodel type specified.');
            throw(err);   
    end

end % end of function


function cmImageGraphPlanted = genBlocks(mDense, vMeanDense, vStdDense, meanSparse, stdSparse)
%
% Generates the block densities for the planted image graph.
%

    cmImageGraphPlanted = cell(1, length(vMeanDense));

    for m = 1 : length(vMeanDense)
    
        mImageGraphPlanted = zeros(size(mDense,1), size(mDense,2)); 

        for c = 1 : size(mDense,2)
            for r = 1 : size(mDense,1)
                if mDense(r,c)
                    density = normrnd(vMeanDense(m), vStdDense(m));
                    if density >= 1.0
                        density = 1.0; 
                    end
                    mImageGraphPlanted(r,c) = density;
                else
                    density = normrnd(meanSparse, stdSparse);
                    if density <= 0.0, density = 0.0; end            
                    mImageGraphPlanted(r,c) = density;
                end
            end
        end
    
        cmImageGraphPlanted{m} = mImageGraphPlanted; % cell of matrixes
    end
    
end



function [cmBackImageMatrix] = genBackgroundImage(cmImageGraphPlanted, vPosNumRow, vPosNumColumn, numRow, numColumn, unitGraphSize)
    %
    % Generate the background image matrix.
    %
    
    posNum = numRow;
    graphSize = max(numRow*unitGraphSize, numColumn*unitGraphSize);
    
    cmBackImageMatrix = cell(1, size(cmImageGraphPlanted,2));
    
    % compute expected number of edges from the cmImageGraphPlaned and the vertex
    % role/position assignments
    for m = 1 : size(cmImageGraphPlanted,2)
        totalEdge = 0;
        for r = 1 : posNum %numRow
            for c = 1 : numColumn
                totalEdge = totalEdge + vPosNumRow(r) * vPosNumColumn(c) * cmImageGraphPlanted{m}(r,c);
            end
        end
        
        % use the total expected edge to obtain the expected edge probability assuming
        % same change of edge
        edgeProb = totalEdge / (graphSize * graphSize);   
        cmBackImageMatrix{m} = ones(numRow, numColumn) * edgeProb;
    end

end


function [cmImageGraph] = genImageGraph(cmImageGraphPlanted, cmImageGraphBackground, vBackgroundProp)
%
% Geneates the image graph from the planted and background image graphs.
%

    assert(size(cmImageGraphPlanted,2) == size(cmImageGraphBackground, 2));

    % each row is a background portion, each column is a generated graph at
    % certain density
    cmImageGraph = cell(length(vBackgroundProp), size(cmImageGraphPlanted,2));
    
    for m = 1 : size(cmImageGraphPlanted,2)
        for b = 1 : length(vBackgroundProp)
            backgroundProp = vBackgroundProp(b);
            cmImageGraph{b,m} = (1-backgroundProp) * cmImageGraphPlanted{m} + backgroundProp * cmImageGraphBackground{m};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for p = 1:size(cmImageGraph{b,m},1)
                for q = 1:size(cmImageGraph{b,m},2)
                    if cmImageGraph{b,m}(p,q) >= 1.0
                        cmImageGraph{b,m}(p,q) = 1.0;
                    end
                    if cmImageGraph{b,m}(p,q) <= 0
                        cmImageGraph{b,m}(p,q) = 0;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


function [vPosNumRow, vPosNumColumn, vVertPos] = genPositionMembership(sPosDist, posNum, numRow, numColumn, unitGraphSize)
%
% Generate the sizes of each position, based on the distribution specified in
% sPosDist.
%
    graphSize = max(numRow*unitGraphSize, numColumn*unitGraphSize);
    switch sPosDist
        case 'uniform'
            vRoleProbRow = zeros(1,numRow) + 1/numRow;
            vRoleProbColumn = zeros(1, numColumn) + 1/numColumn;
        case 'linear'
            vRoleProbRow = [numRow:-1:1]/sum([1:1:numRow]);
            vRoleProbColumn = [numColumn:-1:1]/sum([1:1:numColumn]);
        case 'random'
            % generate a random sizes that sum up to graphsize, then we divide by
            % graphsize
            vRoleProbRow = (randfixedsum(posNum, 1, graphSize, 1, graphSize))' / graphSize;
            vRoleProbColumn = (randfixedsum(numColumn, 1, graphSize, 1, graphSize))' / graphSize;
        otherwise
            disp(usage());
            err = MException('genPositionSize:InvalidPositoinDistribution', 'Invalid role distribution specified.');
            throw(err);
    end
    
    % generate membership of the vertices
    vPosNumRow = mnrnd(numRow*unitGraphSize, vRoleProbRow); 
    vPosNumColumn = mnrnd(numColumn*unitGraphSize, vRoleProbColumn);
    
    % if there are any zero elements, we random sample one of partitions that
    % has more than one element and move an element to each zero element
    % partition
    vZeroSizedRolesRow = find(vPosNumRow == 0);
    for zeroi = 1 : size(vZeroSizedRolesRow,2)
        vNZSizedRolesRow = find(vPosNumRow > 1);
        nzChosen = randsample(vNZSizedRolesRow, 1);   
        vPosNumRow(nzChosen) = vPosNumRow(nzChosen) - 1;
        vPosNumRow(vZeroSizedRolesRow(zeroi)) = vPosNumRow(vZeroSizedRolesRow(zeroi)) + 1;
    end
    vZeroSizedRolesColumn = find(vPosNumColumn == 0);
    for zeroi = 1 : size(vZeroSizedRolesColumn,2)
        vNZSizedRolesColumn = find(vPosNumColumn > 1);
        nzChosen = randsample(vNZSizedRolesColumn, 1);   
        vPosNumColumn(nzChosen) = vPosNumColumn(nzChosen) - 1;
        vPosNumColumn(vZeroSizedRolesColumn(zeroi)) = vPosNumColumn(vZeroSizedRolesColumn(zeroi)) + 1;
    end

    % make sure all roles have at least one element
    assert(isempty(find(vPosNumRow == 0,1)));  
    assert(isempty(find(vPosNumColumn == 0,1)));    
   
    
    % generate vertex to role mapping
    vVertPos = zeros(1, graphSize);
 
    vertIndex = 1;
    for i = 1 : size(vPosNumRow,2)
        % assign vertex, in order, to each role, based on the numbers from
        % vPosNum
        for j = 1 : vPosNumRow(i)
            vVertPos(vertIndex) = i;
            vertIndex = vertIndex + 1;
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cmAdjMat = genGraph(vPosNumRow, vPosNumColumn, cmImageGraph, numRow, numColumn, unitGraphSize)
%
% Generates an adjacency matrix accounting to the postion assignments and given
% image matrix.
%

    % edge value
    edgeVal = 1;

    % adjacency matrix generated
    graphNum = size(cmImageGraph,2);
    posNum = size(cmImageGraph{1},1);
    cmAdjMat = cell(size(cmImageGraph,1), graphNum);
    
    for m = 1 : graphNum
        for b = 1 : size(cmImageGraph,1)
    
            mAdjMat = zeros(numRow*unitGraphSize, numColumn*unitGraphSize); %rectangular
            currRowShift = 1;
            currColShift = 1;

            % from block generate the number of edges (binomial(p,n), p = density of
            % block/probability of an edge, n = number of elements in the block)
            for r = 1 : posNum
                for c = 1 : numColumn
                    edgeNum = binornd(vPosNumRow(r) * vPosNumColumn(c), cmImageGraph{b,m}(r,c)); %a number, like 25752
                    vEdges = randsample(vPosNumRow(r) * vPosNumColumn(c), edgeNum); %(edgeNum * 1) like 36754
                    for e = 1 : size(vEdges,1)
                        row = floor((vEdges(e) - 1) / vPosNumColumn(c));
                        col = floor(mod((vEdges(e)-1), vPosNumColumn(c)));
                        mAdjMat(currRowShift + row, currColShift + col) = edgeVal;
                    end
        
                    currColShift = currColShift + vPosNumColumn(c);
                end
    
                currRowShift = currRowShift + vPosNumRow(r);
                currColShift = 1;
            end
        
            cmAdjMat{b,m} = mAdjMat;
        end
    
    end

end



function [vShuffledVertPos, cmShuffledAdjMat] = shuffle(bShuffle, vVertPos, cmAdjMat, numRow, numColumn, unitGraphSize)
%
% Shuffles the matrices.  The same permutations are applied to all the graphs.
%

    shuffleNum = 50;
    % shuffle adjacency matrix (we need to shuffle the actual partitions also)
    % shuffle x number of times
    
    vShuffledVertPos = vVertPos;
    cmShuffledAdjMat = cmAdjMat;
    graphNum = size(cmAdjMat,2);
    rowSize = numRow*unitGraphSize;
    columnSize = numColumn*unitGraphSize;
    
    if bShuffle
        for s = 1 : shuffleNum
            vShuffleRow = randperm(rowSize);
            vShuffleColumn = randperm(columnSize);
            % shuffle adjacency matrix
            for m = 1 : graphNum
                for b = 1 : size(cmAdjMat,1)
                    cmShuffledAdjMat{b,m} = cmShuffledAdjMat{b,m}(vShuffleRow, vShuffleColumn);
                end
            end
            % shuffle mapping
            vShuffledVertPos = vShuffledVertPos(vShuffleColumn);
        end
    end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sUsage = usage()
%
% Returns the usage string.
%
    sUsage = 'Usage: genBlockmodelWithBackground_rec(sPosDist, posNum, numRow, numColumn, unitGraphSize, bShuffle, vBackgroundProp, sBlockmodelType,varargin)';
end