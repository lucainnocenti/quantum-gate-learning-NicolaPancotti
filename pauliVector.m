function [P] = pauliVector()
    
    sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2);
    P = cell(1,4);
    P{1} = id; P{2} = sx; P{3} = sy; P{4} = sz; 

end

