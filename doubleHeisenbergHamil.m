
function [HHT]=doubleHeisenbergHamil(Jx,Jy,Jz,h,L)
    % L is the length of the chain (L<=19)
    % It will return the full Hamiltonian (sparse matrix)
    % If Jx, Jy, Jz, h are vectors, site dependent 
    % coefficients are used

    d=4;
    D=8;
    sigz=sparse([1,0;0,-1]);
    sigx=sparse([0,1;1,0]);
    sigy=sparse([0,1i;-1i,0]);
    id2=sparse([1,0;0,1]); % eye(2)

    
    xi=kron(sigx,id2);
    ix=kron(id2,sigx.');
    yi=kron(sigy,id2);
    iy=kron(id2,sigy.');
    zi=kron(sigz,id2);
    iz=kron(id2,sigz.');
    II=kron(id2,id2);
    id2 = reshape(id2, [d,1]); 
    
    
    if(L>10)
      warning('Long chain=> slow to compute!!');
       if(L>20)
         error('This L is too big!!')
       end;
    end;

    % Compute "regular" Hamiltonians, operator and vectors

    fprintf('Computing the whole Hamiltonian for L=%d',L);

    Prow = cell(1,D);
    Prow{1} = II; Prow{2} = Jx*xi; Prow{3} = -Jx*ix; Prow{4} = Jy*yi; 
    Prow{5} = -Jy*iy; Prow{6} = Jz*zi; Prow{7} = -Jz*iz; Prow{8} = (h*zi-h*iz);

    Pcol=cell(1,D);
    Pcol{1} = (h*zi-h*iz); Pcol{2} = Jx*xi; Pcol{3} = Jx*ix; Pcol{4} = Jy*yi; 
    Pcol{5} = Jy*iy; Pcol{6} = Jz*zi; Pcol{7} = Jz*iz; Pcol{8} = II;

    B = zeros(d,D,d,D);
    A = zeros(d,1,d,D);
    F = zeros(d,D,d,1);
    
    for k = 1:D
        B(:,1,:,D-k+1) = Prow{k};
        B(:,D-k+1,:,D) = Pcol{k};
        A(:,1,:,k) = Prow{k};
        F(:,k,:,1) = Pcol{k};
    end    
    
    %Special situation: first and last sites are contracted with Id as traced out
    A = permute(A, [1,2,4,3]); A = reshape(A, [d*1*D,d]);
    A = A*id2;
    A = reshape(A, [1,D,d,1]); A = permute(A, [3,1,4,2]);
    
    
    F = permute(F, [1,2,4,3]); F = reshape(F, [d*1*D,d]);
    F = F*id2;
    F = reshape(F, [D,1,d,1]); F = permute(F, [3,1,4,2]);
    
    
    A = reshape(A, [d,D]);
    F = reshape(F, [d,D]);
    
    % 1) on the edges, I contract C1 and CN with their adjoint
    Af = (A')*A;
    Ff = (F')*F;
    
    tmp = B;
    tmp = permute(tmp, [2,1,3,4]); tmp = reshape(tmp, [D,d*D*d]);
    
    edgeL = Af*tmp;
    edgeL = reshape(edgeL, [D*d,1,d,D]);
    
    tmp = B;
    tmp = reshape(tmp, [D*d*d,D]);
    edgeR = tmp*Ff; edgeR = reshape(edgeR,[d,D,d,D]);
    edgeR = permute(edgeR,[1,4,2,3]); edgeR = reshape(edgeR,[d*D,D,d,1]);
    
    
    BT = permute(B,[3,2,1,4]); BT = 1i*BT; BT = reshape(BT,[d,D*d*D]);
    B = permute(B, [1,2,4,3]); B = reshape(B,[D*d*D, d]);
    
    Bf = B*BT;
    Bf = reshape(Bf, [d,D*D,d,D*D]);
    
    edgeLd = BT; edgeLd = reshape(edgeLd, [d,1,D*d,D]); 
    edgeRd = BT; edgeRd = reshape(edgeRd, [d,D,D*d,1]);
   
    edgeLd = permute(edgeLd, [1,2,4,3]); edgeLd = reshape(edgeLd, [d*D,d*D]);
    edgeRd = permute(edgeRd, [1,2,4,3]); edgeRd = reshape(edgeRd, [d*D,d*D]);
    
    edgeR = permute(edgeR,[4,2,3,1]); edgeR = reshape(edgeR,[D*d*1,d*D]);
    edgeL = permute(edgeL,[4,2,3,1]); edgeL = reshape(edgeL,[D*d*1,d*D]);
    
    
    
    edgeA = edgeL*edgeLd;
    edgeF = edgeR*edgeRd;
    
    edgeA = reshape(edgeA, [d,1,d,D*D]);
    edgeF = reshape(edgeF, [d,D*D,d,1]);
    
    HHT = cell(1,L);
    
    HHT{1} = edgeA; HHT{L} = edgeF;
    
    for j = 2:(L-1) 
        HHT{j} = Bf;
    end
    
    

end