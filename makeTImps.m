function [mps] = makeTImps(N, A, F, B)
    % This function creates a translationally invariant mps. N is the
    % number of sites, A is the initial site tensor, F is the final and B 
    %is the inner.
    %Usually, for mps of spin systems, A,F are 2xD matrices where 2 is the 
    %physical dymension and D the bond dimension. While B is a 2xDxD tensor
    %with the same meaning.
    
    mps = cell(1,N);
    
    mps{1} = A/sqrt(4); mps{N} = F/sqrt(4);
    
    for j = 2:(N-1)
        
        mps{j} = B/sqrt(4);
        
    end

end

