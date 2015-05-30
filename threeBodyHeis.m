function [H3] = threeBodyHeis (Jx,Jy,Jz,h)

    sigz=sparse([1,0;0,-1]);
    sigx=sparse([0,1;1,0]);
    sigy=sparse([0,-1i;1i,0]);
    id2=sparse([1,0;0,1]); % eye(2)

    xx=kron(sigx,sigx);
    yy=kron(sigy,sigy);
    zz=kron(sigz,sigz);
    zi=kron(sigz,id2);
    iz=kron(id2,sigz);
    
    
    H3 = Jx*kron(xx,id2) + Jy*kron(yy,id2) + Jz*kron(zz,id2);
    H3 = H3 + Jx*kron(id2,xx) + Jy*kron(id2,yy) + Jz*kron(id2,zz);
    H3 = H3 + h*kron(zi,id2) + h*kron(iz,id2) + h*kron(id2,iz);