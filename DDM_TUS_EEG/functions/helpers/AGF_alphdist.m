function [ dist ] = AGF_alphdist( A, B )
%Little function that determines the distance of two letters in the
%alphabet from each other.
    ABC='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    Ap=find(ABC== upper(A));
    Bp=find(ABC== upper(B));
    dist=abs(Ap-Bp);
end

