function oper=expandOper(Ops,d,N)
  % Take array of N operators and expand to d^N x d^N single operator
  tmp=1;
  for k=1:N
    disp(sprintf('Expanding pos %d',k))
    %disp(Ops{k})
    [dphys1,dphys2,Dl]=size(tmp);
    [d1,Dlp,d2,Dr]=size(Ops{k});
    if(Dlp~=Dl)
    disp(sprintf('Failed at Step %d, tmp(%d,%d,%d) Multiplying next operator(%d,%d,%d,%d)',k,...
                 size(tmp,1),size(tmp,2),size(tmp,3),d1,Dlp,d2,Dr))
        error(sprintf('Dimensions do not agree in expandOper! %d != %d',Dlp,Dl));
    end;
    disp(sprintf('Step %d, size=%d,%d,%d Multiplying next operator',k,...
                 size(tmp,1),size(tmp,2),size(tmp,3)))
    tmp=reshape(tmp,dphys1*dphys2,Dl)*...
        reshape(permute(Ops{k},[2 1 3 4]),Dlp,d1*d2*Dr);
    disp(sprintf('Reshaping to %d x %d x %d',dphys1*d1,dphys2*d2,Dr))
    tmp=reshape(permute(reshape(tmp,dphys1,dphys2,d1,d2,Dr),[1 3 2 4 5]),...
                dphys1*d1,dphys2*d2,Dr);
  end;
  oper=tmp;
  disp('End of expandOper')
  end