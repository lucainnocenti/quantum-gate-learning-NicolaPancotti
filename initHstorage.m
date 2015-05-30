function [Hstorage]=initHstorage(mps,hset,d)

    %%%%%%%%%%%%%%%%%%%%%%
    %%First create M by N empty cells
    %%M is the number of interactions
    %%N is the number of spins/particles
    %%For example, an Heisenberg-like hamiltonian with no magnetic field
    %%has two body diagonal interactions (xx,yy,zz). In this case M=3(N-1)
    %%%%%%%%%%%%%%%%%%%%%
    [M,N]=size(hset);
    Hstorage=cell(M,N+1);
    %%We impose the boundary conditions 
    for m=1:M, Hstorage{m,1}=1; Hstorage{m,N+1}=1; end
    %%We contract each cell with the corresponding hamiltonian item. The
    %%j-th element of Hstorage() contains the contraction of all the
    %%tensors from N to j with the respective hamiltonian element 
    for j=N:-1:2
        for m=1:M
            h=reshape(hset{m,j},[1,1,d,d]);
            Hstorage{m,j}=updateCright(Hstorage{m,j+1},mps{j},h,mps{j});
        end
    end
end

