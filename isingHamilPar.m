function [H,H1,H2]=isingHamilPar(L,J,g,h)
%> Return the following Hamiltonian
%> H=Sum_i -J[sigmaz_i sigmaz_i+1 + g/2 sigmax_i + g/2 sigmax_i+1
%>                                + h/2 sigmaz_i + h/2 sigmaz_i+1]
% Also one-body (H1) and 2-body (H2) pieces, to be able to compare to
% Trotter expansions
 
sigz=sparse([1,0;0,-1]);
sigx=sparse([0,1;1,0]);
id2=sparse([1:2],[1:2],1); % eye(2)

zz=kron(sigz,sigz);
zi=kron(sigz,id2);
xi=kron(sigx,id2);

d=2;

if(L>10)
  warning('Long chain=> slow to compute!!');
   if(L>20)
     error('This L is too big!!')
   end;
end;

% Compute "regular" Hamiltonian
disp(sprintf('Computing the whole Ising Hamiltonian for L=%d',L))

H=sparse(d^L,d^L);H1=H;H2=H;
for k=1:L-1
    %locT=J*zz+g*xi+h*zi;    
    idL=speye(d^(k-1));
    idR=speye(d^(L-k-1));
    
    auxH2=J*kron(idL,kron(zz,idR));
    auxH1=kron(idL,kron(g*xi+h*zi,idR));
    H=H+auxH1+auxH2;
    H1=H1+auxH1;
    H2=H2+auxH2;
end;

% Plus the last single term
lastT=g*kron(speye(d^(L-1)),sigx)+h*kron(speye(d^(L-1)),sigz);
H=H+lastT;
H1=H1+lastT;
