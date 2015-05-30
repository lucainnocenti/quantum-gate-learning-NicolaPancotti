function [hset] = HeisenbergHamiltonian(N)
    % Heisenberg Hamiltonian
    M=3*(N-1);
    hset=cell(M,N);
    sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2);
    for m=1:M, for j=1:N, hset{m,j}=id; end; end
    for j=1:(N-1)
        hset{3*(j-1)+1,j}=sx; hset{3*(j-1)+1,j+1}=sx;
        hset{3*(j-1)+2,j}=sy; hset{3*(j-1)+2,j+1}=sy;
        hset{3*(j-1)+3,j}=sz; hset{3*(j-1)+3,j+1}=sz;
    end

end

