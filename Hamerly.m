set=input('Odaberi set podataka: \n 1: uniform random \n2:uniform random \n3: r15 - \n 4: pathbased - \n 5: covtype\n 6: kddcup\n');
switch set
    case 1
        A=rand(50000,2); 
    case 2
        A=rand(50000,8);
    case 3
        A=load('r15.txt');
        A=A(:,1:2); %treba ispast n=600, pathbased-> 300
    case 4
        A=load('pathbased.txt');
        A=A(:,1:2);
    case 5
        A=load('covertype.mat');
        A=A.B;
    case 6
        filename = 'cup98.txt';
        delimiterIn = ' ';
        headerlinesIn = 1;
        A = importdata(filename,delimiterIn,headerlinesIn);
        A=A.data;
    otherwise
        display('Izabrat cu sama');
        A=load('pathbased.txt');
        A=A(:,1:2);
        
end
k=input('Unesi broj klastera:\n');
n=size(A,1);
dim=size(A,2);
centri=zeros(k,dim);
flag=input('Odaberi jedno od sljedeceg:\n 1: Koristi gotove centre\n 2: Isprobaj nove (pa mozes onda i za drukciji broj klastera) - random metoda\n 3: further-first inicijalizacija\n' );
if flag==1
    filename='ime datoteke';
    filename=input('Upisi ime datoteke s centrima koju trebas, npr. centers3-3.mat\n \n', 's');
    centri=load(filename);
    centri=centri.centri; 
elseif flag==2
    ind=randperm(n);
    ind=ind(1:k);
    for j=1:k
        centri(j,:)=A(ind(j),:);
    end
else 
   for i=1:n
       centri(1,:)=centri(1,:)+A(i,:);
   end
   centri(1,:)=centri(1,:)/n;
   for j=2:k
       
       maksi=zeros(1,n);
       for i=1:n
           dist=zeros(1,j-1);
           for l=1:j-1
               dist(l)=norm(A(i,:)-centri(l,:));
           end
           maksi(i)=min(dist);
       end
       [maxi,maxarg]=max(maksi);
       centri(j,:)=A(maxarg,:);
   end
    
end
dim=size(A,2);
k=size(centri,1);
n=size(A,1);

indeksi=zeros(1,n);
q=zeros(1,k);
p=zeros(1,k);
s=zeros(1,k);
u=inf*ones(1,n);
l=zeros(1,n);
br_praznih=0;
cluster_sum=zeros(k,dim);
br_iter=1;
promjena=1;
%inicijalizacija
%postavljanje q i cluster_sum na nule
for i=1:n
   %PointAllCtrs
    dist=zeros(1,k);
    for j=1:k
        dist(j)=norm(A(i,:)-centri(j,:));
    end
    [mini, minarg]=min(dist);
    indeksi(i)=minarg; %minarg je indeks klustera za koji se postiže najmanja udaljenost
    u(i)=mini; %udaljenost do najblizeg centra, u(i)=dist(indeksi(i));
    dist(indeksi(i))=inf; %ovdje bi trebala razmislit u slucaju da je k=1
    [mini, minarg]=min(dist);
    l(i)=mini; %udaljenost do drugog najblizeg centra
    q(indeksi(i))=q(indeksi(i))+1;
    cluster_sum(indeksi(i),:)=cluster_sum(indeksi(i),:)+A(i,:);
end


%pomicem centre
temp=centri;
    for j=1:size(centri,1)

        if ~q(j) 
            br_praznih=br_praznih+1;
            [maksi,maxarg]=max(q); %trazim najveci kluster i njemu vadim tocku
            trazi=find(indeksi==maxarg);
            ind=trazi(1); %indeks prvog po redu elementa u najvecem klusteru
            centri(j,:)=A(ind,:);
            indeksi(i)=j;
            cluster_sum(j,:)=A(ind,:);
            q(j)=1;
            cluster_sum(maxarg,:)=cluster_sum(maxarg,:)-A(ind,:);
            q(maxarg)=q(maxarg)-1;
            if maxarg<j %ako smo u prethodnoj iteraciji ovo racunali, moramo ponovo
                centri(maxarg,:)=cluster_sum(maxarg,:)/q(maxarg);
                p(j)=norm(temp(maxarg,:)-centri(maxarg,:));
            end
            
   
        else
            centri(j,:)=cluster_sum(j,:)/q(j); %OVOOOOO cluster_sum(j)
            p(j)=norm(temp(j,:)-centri(j,:));
        end
    end
    

  [maksi, maxarg]=max(p);
    r=maxarg;
    temp=p;
    temp(maxarg)=-inf;
    [maksi, maxarg]=max(temp);
    r2=maxarg;

    for i=1:size(u,2)
        u(i)=u(i)+p(indeksi(i));
        if r==indeksi(i)
            l(i)=l(i)-p(r2);
        else
            l(i)=l(i)-p(r);
        end

    end
    
    
    
br_ent=0;
tic
while promjena
    promjena=0;
    br_iter=br_iter+1;
    for j=1:size(centri,1)
        dist=zeros(1,size(centri,1));
        for r=1:size(centri,1)
            dist(r)=norm(centri(j,:)-centri(r,:));
        end
        [mini, minarg]=min(dist);
        dist(minarg)=inf; %jer je minarg=j pa je dist(minarg)=0, a mi zelimo najblizi centar ovome
        [mini, minarg]=min(dist);
        s(j)=mini;
    end
    
    for i=1:size(A,1)
        m=max(s(indeksi(i))/2, l(i));
        if u(i)>m
            u(i)=norm(A(i,:)-centri(indeksi(i),:)); 
            if u(i)>m
                br_ent=br_ent+1; %broj ulaska u najdublju petlju
                indeksi2=indeksi(i);
                %PointAllCtrs
                dist=zeros(1,size(centri,1));
                l2=inf;
                for j=1:size(centri,1)
                    dist(j)=norm(A(i,:)-centri(j,:));
                    if dist(j)<u(i)
                        l2=u(i);
                        u(i)=dist(j);
                        indeksi(i)=j;
                    elseif dist(j)<l2
                        l2=dist(j);
                    end
                end
                l(i)=l2;
                if indeksi2~=indeksi(i)
                    %trazimo broj tocaka u klusteru a'=indeksi2 -->
                    %prebrojimo koliko se puta a' pojavljuje u nizu a
                    q(indeksi2)=q(indeksi2)-1;
                    q(indeksi(i))=q(indeksi(i))+1;
                    cluster_sum(indeksi2,:)=cluster_sum(indeksi2,:)-A(i,:);
                    cluster_sum(indeksi(i),:)=cluster_sum(indeksi(i),:)+A(i,:);
                    
                end
            else continue;
            end
            
        else continue;
        end;
    end
    
    %MoveCenters
    temp=centri;
    for j=1:size(centri,1)

        if ~q(j) 
            br_praznih=br_praznih+1;
            [maksi,maxarg]=max(q); %trazim najveci kluster i njemu vadim tocku
            trazi=find(indeksi==maxarg);
            ind=trazi(1); %indeks prvog po redu elementa u najvecem klusteru
            centri(j,:)=A(ind,:);
            indeksi(i)=j;
            cluster_sum(j,:)=A(ind,:);
            q(j)=1;
            cluster_sum(maxarg,:)=cluster_sum(maxarg,:)-A(ind,:);
            q(maxarg)=q(maxarg)-1;
            if maxarg<j %ako smo u prethodnoj iteraciji ovo racunali, moramo ponovo
                centri(maxarg,:)=cluster_sum(maxarg,:)/q(maxarg);
                p(j)=norm(temp(maxarg,:)-centri(maxarg,:));
            end
            
   
        else
            centri(j,:)=cluster_sum(j,:)/q(j);
            p(j)=norm(temp(j,:)-centri(j,:));
        end
    end
    [maksi, maxarg]=max(p);
    if maksi>0 
        promjena=1;
    else
        break;
    end
    %UpdateBounds
    [maksi, maxarg]=max(p);
    r=maxarg;
    temp=p;
    temp(maxarg)=-inf;
    [maksi, maxarg]=max(temp);
    r2=maxarg;

    for i=1:size(u,2)
        u(i)=u(i)+p(indeksi(i));
        if r==indeksi(i)
            l(i)=l(i)-p(r2);
        else
            l(i)=l(i)-p(r);
        end

    end
end


toc

display(br_iter);

rez=br_ent/(br_iter*n);
display(rez);

x=input('Zelite li plot?\n 0: ne\n 1: da\n');
if x
 for i=1:n
     plot(A(i,1), A(i,2), '.b');
     if(i==1) hold on;
     end
 end
 
 for i=1:k
     plot(centri(i,1), centri(i,2), '+r');
end
hold off
end


