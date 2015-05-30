function [H,He,Ho,op0]=HeisenbergHamil(Jx,Jy,Jz,h,L)
% L is the length of the chain (L<=19)
% It will return the full Hamiltonian (sparse matrix)
% If Jx, Jy, Jz, h are vectors, site dependent 
% coefficients are used
 
d=2;
sigz=sparse([1,0;0,-1]);
sigx=sparse([0,1;1,0]);
sigy=sparse([0,-i;i,0]);
id2=sparse([1:2],[1:2],1); % eye(2)

xx=kron(sigx,sigx);
yy=kron(sigy,sigy);
zz=kron(sigz,sigz);
zi=kron(sigz,id2);
iz=kron(id2,sigz);

if(L>10)
  warning('Long chain=> slow to compute!!');
   if(L>20)
     error('This L is too big!!')
   end;
end;

% Coefficients
if(length(Jx)<L-1)
    % repeat the first value until the end
    Jx(length(Jx)+1:L-1)=Jx(1);
end;
if(length(Jy)<L-1)
    % repeat the first value until the end
    Jy(length(Jy)+1:L-1)=Jy(1);
end;
if(length(Jz)<L-1)
    % repeat the first value until the end
    Jz(length(Jz)+1:L-1)=Jz(1);
end;
if(length(h)<L)
    % repeat the first value until the end
    h(length(h)+1:L)=h(1);
end;


% Compute "regular" Hamiltonians, operator and vectors

disp(sprintf('Computing the whole Hamiltonian for L=%d',L))

H=sparse(d^L,d^L);
He=sparse(d^L,d^L);Ho=sparse(d^L,d^L);
for k=1:L-1
    locT=Jx(k)*xx+Jy(k)*yy+Jz(k)*zz+h(k)*zi;    
    locH12=Jx(k)*xx+Jy(k)*yy+Jz(k)*zz+.5*h(k)*zi+.5*h(k+1)*iz;
    if(k==1) locH12=locH12+.5*h(k)*zi; end;
    if(k==L-1) locH12=locH12+.5*h(k+1)*iz; end;
    idL=speye(d^(k-1));
    idR=speye(d^(L-k-1));
    auxH=kron(idL,kron(locT,idR));
    H=H+auxH;
    if(mod(k-1,2)==0)
        He=He+kron(idL,kron(locH12,idR));
    else
        Ho=Ho+kron(idL,kron(locH12,idR));
    end;
end;

% Plus the last single term
lastT=h(L)*kron(speye(d^(L-1)),sigz);
H=H+lastT;