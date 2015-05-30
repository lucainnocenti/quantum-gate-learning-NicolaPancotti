function [MPO, MPS] = permuteOpers(mpo, mps)

    if size(mps,2) ~= size(mpo,2)
        error('mpo & mps don t have the same length')
    end;
    
    MPO = cell(1,size(mpo,2));
    MPS = cell(1,size(mpo,2));
        
    for i = 1:size(mpo,2)
        
        MPO{i} = permute(mpo{i},[2,4,1,3]);
        MPS{i} = permute(mps{i},[2,3,1]);
        
    end