tic
[A, map]=imread('mala_loptica.jpg');
if ~isempty(map)
    Im=ind2rgb(A,map);
end

M=im2double(A); %A sadrzi podatke tipa uint8 pa konvertiramo u double

[r,c, ~]=size(M); 
l=r*c*3; 

Dataset=reshape(M,[r*c,3]); %svaki redak predstavlja piksel, stupci rgb vrijednosti

n=r*c; %broj podataka u skupu
%sad radimo matricu afiniteta
W=zeros(n,n);
for i=1:n-1
    for j=i+1:n
%        piksel koji je u originalnoj matrici bio na poziciji (k, l), sad je
%        na poziciji r*(l-1)+k=i --> k=mod(i, r); l=ceil(i/r);
%         a=mod(i,r); b=ceil(i/r);
%         c=mod(j,r); d=ceil(j/r);
%         if a==0 
%            a=r;
%         end
%         if c==0 
%             c=r;
%         end
        %(a,b) koordinate piksela u i-tom retku, (c,d) piksela u j-tom
        %retku
        %dist=sqrt((a-c)^2+(b-d)^2);
        boje=10*norm(Dataset(i,:)-Dataset(j,:)); 
        W(i,j)= exp(-boje^2);
    end
end

W=W+W';
 nColors=3;
X=spectral_clustering(W,nColors, eps);

for i=1:n
     a=mod(i,r);
     if a==0
         a=r;
     end
     %a je redak
     b=ceil(i/r); %stupac
     
     if X(i, 1) == 1
        M(a,b,1)=1;
        M(a,b,2)=0;
        M(a,b,3)=0;
        
     elseif X(i,2)==1
         M(a,b,1)=0;
        M(a,b,2)=1;
        M(a,b,3)=0;
     elseif X(i, 3) == 1
         M(a,b,1)=0;
         M(a,b,2)=0;
         M(a,b,3)=1;
         
     end
     

end

figure(2), imshow(M,map)
figure(1), imshow(A,map)
toc