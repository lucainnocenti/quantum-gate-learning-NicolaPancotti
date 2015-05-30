function [e,n]=expectationvalue2(mps,hset)
    
    [M,N]=size(hset);
    d=size(mps{1},3);
    % expectation value
    e=1;
    for j=N:-1:1
        h=hset{j};
        e=updateCright(e,mps{j},h,mps{j});
    end
    
   
    % norm
    n=1;
    X=eye(d); X=reshape(X,[1,1,d,d]);
    for j=N:-1:1
        n=updateCright(n,mps{j},X,mps{j});
    end
    e=e/n;
end
