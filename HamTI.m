function [hset,B] = HamTI(N)
    
    P = pauliVector();
    Prow = [P,0];
    temp1 = Prow(1);
    temp5 = Prow(5);
    Pcol = Prow;
    Pcol(1) = temp5;
    Pcol(5) = temp1;
    
    B = zeros(5,5,2,2);
    for i = 1:size(Prow,2)
        B(1,i,:,:) = Prow{i};
    end
    for i = 1:size(Pcol,2)
        B(i,5,:,:) = Pcol{i};
    end
    A = zeros(1,5,2,2);
    for i = 1:size(Prow,2)
        A(1,i,:,:) = Prow{i};
    end
    F = zeros(5,1,2,2);
    for i = 1:size(Pcol,2)
        F(i,1,:,:) = Pcol{i};
    end
        
        
    hset = cell(1,N);
    
    hset{1} = A;
    hset{N} = F;
    
    for j = 2:(N-1)
        hset{j} = B;
    end
    
end

