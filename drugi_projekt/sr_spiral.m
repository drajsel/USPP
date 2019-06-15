A=load('spiral.txt');
A=A(:,1:2);
k=3;

% redovi od A su tocke

%konsturiranje matrice W
n=size(A,1);
W=zeros(n,n);
brojac=0;
sum=0;
for i=1:n
   for j=i+1:n
       sum=sum+norm(A(i,:)-A(j,:));
       brojac=brojac+1;
   end
end

sigma=sum/brojac;
       
for i=1:n
    for j=1:n
        W(i,j)=exp(-norm(A(i,:)-A(j,:))^2/sigma);
    end
end

jedinice=ones(n,1);
D=diag(W*jedinice);
XT=zeros(n,k);
X2=zeros(n,k);
%korijen=inv(sqrt(D)), tj. D^(-1/2)
%korijen=zeros(n,n); 
%korijen=D^(-1/2);

korijen=D;
for i=1:n
    korijen(i,i)=1/sqrt(korijen(i,i));
end
H=korijen*W*korijen; %na nju primjenjujemo eig u 2. koraku algoritma

[V, E]=eigs(H,k); %E je S iz algoritma, H mora bit simetricna

% [Y,E]=eig(H);
% for i=1:k
%     V(:,i)=Y(:,n-k+i); % V sadrzi svojstvene vektore k najvecih svojstvenih vrijednosti
% end

Z=korijen*V; %nxk matrica
%normaliziramo Z, normiramo retke od Z i to nam je X tilda
for i=1:n
    XT(i,:)=Z(i,:)/norm(Z(i,:));
end
%XT=diag(1./sqrt(diag(Z*Z')))*Z;
%cetvrti korak algoritma, inicijalizacija
R=zeros(k,k);
i=randi(n,1);
R(:,1)=XT(i,:)';
c=zeros(n,1);
for l=2:k
    c=c+abs(XT*R(:,l-1));
    [mini, minarg]=min(c);
    R(:,l)=XT(minarg,:)';
end

conv=0;
X=zeros(n,k);
br_iter=0;
max_it=10000;
iter=zeros(1,max_it);
while br_iter<max_it
    br_iter=br_iter+1;
    %diskretiziramo rjesenje (6. korak)
    X2=XT*R; %mi smo imali XT zapravo X tilda sa zvjezdicom,
    %a X2 je Xtilda iz koraka 6.
    for i=1:n
        [maxi, maxarg]=max(X2(i,:));
        X(i,:)=zeros(1,k);
        X(i,maxarg)=1;
    end
   %trazimo ortonormiranu R (7. korak)
   %napravimo SVD od X'XT
   [P,S,Q]=svd(X'*XT);
   
   phi=trace(S);
   iter(br_iter+1)=abs(phi-conv);
   if abs(phi-conv)<eps
       break;
   end
    conv=phi;
    R=Q*P';
end
  
centers=zeros(k,2);
zastava=0;
for l=1:k
    brojac=0;
    for i=1:n
        if X(i,l)
            centers(l,:)=centers(l,:)+A(i,:);
            brojac=brojac+1;
            if ~zastava
                zastava=1;
                hold on;
            end
            if l==1
                plot(A(i,1), A(i,2), '*b', 'Markersize', 20);
            elseif l==2
                plot(A(i,1), A(i,2), '*r', 'Markersize', 20);
            else 
                plot(A(i,1), A(i,2), '*g', 'Markersize', 20);
            
            end
        end
    end
    centers(l,:)=centers(l,:)/brojac;
    
end
hold off;


