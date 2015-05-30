function [ y ] = commutator(x)
    
    H = HeisenbergHamil(1,1,1,3,6);
    
    Hi = kron(eye(2^6),H) - kron(H,eye(2^6));
    HiiH = Hi*Hi';
    
    y = HiiH*x;
    
     

end

