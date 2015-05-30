function [X,numindX]=contracttensors(X,numindX,indX,Y,numindY,indY)
   
    %This function contract the indeces indX and inY of tensor X and Y
    %respectively
    
    %%%%%%%%%%%%%%%%%%
    % sizeX/sizeY contain the indeces of X/Y
    % indX/indY contain the indeces to contract
    % sizeX/sizeY contain the size in each index in indX/indY 
    % indXl/indYr contain the indeces to not contract
    % sizeXl/sizeYr contain the size in each index in indXl/indYr
    %%%%%%%%%%%%%%%%%%

    Xsize=ones(1,numindX); Xsize(1:length(size(X)))=size(X);
    Ysize=ones(1,numindY); Ysize(1:length(size(Y)))=size(Y);
    
    indXl=1:numindX; indXl(indX)=[];
    indYr=1:numindY; indYr(indY)=[];
    
    sizeXl=Xsize(indXl);
    sizeX=Xsize(indX);
    sizeYr=Ysize(indYr);
    sizeY=Ysize(indY);
    
    if prod(sizeX)~=prod(sizeY)
        error('indX and indY are not of same dimension.');
    end
    
    if isempty(indYr)
        if isempty(indXl)
            
            % In case we contract all the indeces
            
            X=permute(X,[indX]);
            X=reshape(X,[1,prod(sizeX)]);
            
            Y=permute(Y,[indY]);
            Y=reshape(Y,[prod(sizeY),1]);
            
            X=X*Y;
            Xsize = 1;
            
        return;
        
        else
            
            %in case we contract all the indeces of the tensor Y
            
            X=permute(X,[indXl,indX]);
            X=reshape(X,[prod(sizeXl),prod(sizeX)]);
            Y=permute(Y,[indY]);
            Y=reshape(Y,[prod(sizeY),1]);
            X=X*Y;
            Xsize=Xsize(indXl);
            X=reshape(X,[Xsize,1]);
            
            return
        end
    end
    
    X=permute(X,[indXl,indX]); % I put on the left the indeces I don't contract  
    X=reshape(X,[prod(sizeXl),prod(sizeX)]); % I transform the tensor into a matrix
    
    Y=permute(Y,[indY,indYr]); % I put on the right the indeces I don't contract
    Y=reshape(Y,[prod(sizeY),prod(sizeYr)]);
    
    X=X*Y; % Perform the multiplication between the X right and Y left
    Xsize=[Xsize(indXl),Ysize(indYr)]; %define the size of the remaining object
    
    numindX=length(Xsize); % ...the number of its indeces
    X=reshape(X,[Xsize,1]);% why [Xsize,1] ???
    
    
end
