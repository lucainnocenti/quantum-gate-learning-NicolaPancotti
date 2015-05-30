function [MPOS] = reshapeMps(mps)

    MPOS = cell(1,size(mps,2));
        
    for i = 1:size(mps,2)
        
        MPOS{i} = reshape(mps{i},[2,2,size(mps{i},2),size(mps{i},3)]);
        MPOS{i} = permute(MPOS{i},[1,3,2,4]);
        
    end