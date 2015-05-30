function [hset,B] = HeisHamTI(N,h)
        
    sz=sparse([1,0;0,-1]);
    sx=sparse([0,1;1,0]); 
    sy=sparse([0,-1i;1i,0]);    
    id=sparse([1,0;0,1]);
    Prow = cell(1,5);
    Prow{1} = id; Prow{2} = sx; Prow{3} = sy; Prow{4} = sz; Prow{5} = h*sz; 
    
    
    Pcol = Prow;
    temp2 = Prow(5);
    temp1 = Prow(1);
    Pcol(1) = temp2;
    Pcol(5) = temp1;
    
    A = zeros(1,5,2,2);
    for j = 1:size(Prow,2)
        A(1,j,:,:) = Prow{j};
    end
    F = zeros(5,1,2,2);
    for j = 1:size(Pcol,2)
        F(j,1,:,:) = Pcol{j};
    end
    
    t = cell(1); t{1} = [0,0;0,0]; 
    Pcol(1) = t(1);
    
    B = zeros(5,5,2,2);
    
    for j = 1:size(Pcol,2)
        B(j,5,:,:) = Pcol{j};
    end
    for j = 1:size(Prow,2)
        B(1,j,:,:) = Prow{j};
    end
    
    
        
        
    hset = cell(1,N);
    
    hset{1} = A;
    hset{N} = F;
    
    for j = 2:(N-1)
        hset{j} = B;
    end
    
end

