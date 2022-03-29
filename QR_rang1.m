function [Q1,R1,A1]=QR_rang1(A,u,v)
%% Actualizarea factorizarii QR cand ea sufera o modificare de rang 1
%  INPUTS
%   A   --  matrice aleatoare de dimensiune (m,n),
%   u   --  vector aleator de dimensiune (m,1),
%   v   --  vector aleator de dimensiune (n,1).
%
%  OUTPUTS
%   Q1   --  matrice ortogonala de dimensiune (m,n),
%   R1   --  matrice superior triunghiulara de dimensiune (m,n),
%   A1   --  matricea A dupa ce a suferit modificarea de rang 1 (A1=A+u*v').
%
%% SOLUTION START %%
[Q,R]=qr(A);
[m,n]=size(A);   % dimensiunile lui A
A1=A+u*v';     %calculez A1, dupa ce A a suferit modificarea
w=Q'*u;    % calculez w din ecuatia lui A1

for i=m-1:-1:1    % pentru a calcula cele m-1 rotatii Jk
    
    [G,y]=planerot([w(i) w(i+1)]');  %calculez o rotatie astfel incat sa introduc un 0 in w
    %extrag din G, c si s
    c=G(1,1);
    s=-G(1,2); % pentru ca lucrez cu rotatiile de tip [c -s; s c] a trebuit sa inmultesc cu -1
    
    w(i+1)=0;
    w(i)=y(1);   % ia valoarea normei lui w(i),w(i+1)
    
    for j=n:-1:1 % parcurg coloanele lui R si aplic rotatia calculata anterior
        alfa=R(i,j);
        R(i,j) = c*R(i,j) - s*R(i+1,j);
        R(i+1,j) = s*alfa + c*R(i+1,j);
    end
    for j=m:-1:1   % parcurg liniile lui Q si aplic rotatia calculata anterior
        alfa=Q(j,i);
        Q(j,i) = c*Q(j,i) - s*Q(j,i+1);
        Q(j,i+1) = s*alfa + c*Q(j,i+1);
    end
end
% dupa iesirea din bucla R va avea forma sup. Hessenberg, iar w va avea un
% singur element nenul pe prima linie
R(1,:)=R(1,:)+w(1)*v';

for k=1:n-1   %pentru a calcula cele n-1 rotatii Gk
    [G,y]=planerot([R(k,k) R(k+1,k)]');  %calculez o rotatie astfel incat sa introduc sub diagolona principala a lui R un 0
    c=G(1,1);
    s=-G(1,2);
    
    for j=1:n     %parcurg coloanele lui R si aplic rotatia calculata anterior
        alfa=R(k,j);
        R(k,j)=c*R(k,j)-s*R(k+1,j);
        R(k+1,j)=s*alfa+c*R(k+1,j);
    end
    for j=1:m   %parcurg liniile lui Q si aplic rotatia calculata anterior
        alfa=Q(j,k);
        Q(j,k)=c*Q(j,k)-s*Q(j,k+1);
        Q(j,k+1)=s*alfa+c*Q(j,k+1);
    end
end

Q1=Q;
R1=R;
%% SOLUTION END %%
end