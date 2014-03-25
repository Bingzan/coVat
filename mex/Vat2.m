function [mRearrangedDis, vPermVerts, mMst] = Vat2(mDis)
%
% Wrapper for c implementation of VAT.
%
% @author: Jeffrey Chan, 2014
%
    [vPermVerts, mMst] = cvat(mDis);
    mMst = mMst + mMst';
    mRearrangedDis = mDis(vPermVerts, vPermVerts);

end % end of function