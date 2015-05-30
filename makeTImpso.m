function [mpo] = makeTImpso(N, A, F, B)
    % This function creates a translationally invariant mpo. N is the
    % number of sites, A is the initial site tensor, F is the final and B 
    %is the inner.
    %Usually, for mpo of spin systems, A,F are 2x2xD matrices where 2 is the 
    %physical dymension and D the bond dimension. While B is a 2x2xDxD tensor
    %with the same meaning.
    
    mpo = cell(1,N);
    
    mpo{1} = A; mpo{N} = F;
    
    for j = 2:(N-1)
        
        mpo{j} = B;
        
    end

end
