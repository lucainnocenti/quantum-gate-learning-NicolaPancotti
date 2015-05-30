function [mpsOr] = ordering(mps)

    L = size(mps,2);
    
    mpsOr = cell(1,L);
    
    for i= 1:L
        
        mps{i} = permute(mps{i}, [1,3,2]);
        mpsOr{L+1-i} = mps{i};
        
    end
    