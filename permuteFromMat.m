function [MPO] = permuteFromMat(mpo)

    
    MPO = cell(1,size(mpo,2));
        
    for i = 1:size(mpo,2)
        
        MPO{i} = permute(mpo{i},[3,1,4,2]);
      
        
    end