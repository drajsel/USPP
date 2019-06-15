
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
novi_centri=zeros(k,dim);
l=zeros(n,k); % matrica donjih ograda na udaljenost tocke do centara
u=inf*ones(1,n); % niz gornjih ograda na udaljenost tocke od centara
s=zeros(1,k); % s(c) je polovica udaljenosti centra od drugog najblizeg centra
dist_center=zeros(k,k);
%dist_center je matrica udaljenosti za centre
% for i=1:k
%     for j=1:k 
%         dist_center(i,j)=norm(centri(i,:)-centri(j,:));
%     end
% end
for i=1:n 
    indeksi(i)=randi(k,1);
end  
%u prvoj iteraciji svakoj toèki dodijelimo najbliži centar
for i=1:n
    dist=zeros(1,k);
    for j=1:k
            dist(j)=norm(centri(j,:)-A(i,:));
    end
    [mini, minarg]=min(dist);
    indeksi(i)=minarg;
end
brojac=zeros(1,k);
for j=1:k
    centri(j,:)=zeros(1,dim);
    brojac(j)=0;
    
    for i=1:n
        if indeksi(i)==j
            centri(j,:)=centri(j,:)+A(i,:);
            brojac(j)=brojac(j)+1;
        end
    end
    centri(j,:)=centri(j,:)/brojac(j);
end
promjena=1;
br_iter=1;
br_praznih=0;
br_ent=0;
tic
while(promjena)
    promjena=0;
    br_iter=br_iter+1;
        for i=1:k
            for j=1:k
            dist_center(i,j)=norm(centri(i,:)-centri(j,:));
            end
        end
 
     %sad racunamo najmanju udaljenost dva razlicita centra pa podijelimo s
     %2
    
    for i=1:k
        mini=Inf;
        for j=1:k
            if(j~=i && dist_center(i, j)<mini) mini=dist_center(i, j);
            end
        end
        s(i)=0.5*mini;
    end

    for i=1:n
        if(u(i)>s(indeksi(i)))
            br_ent=br_ent+1;
            for j=1:k
                if(j~=indeksi(i) && u(i)>l(i,j) && u(i)>0.5*dist_center(indeksi(i), j))
                dist=norm(A(i,:)-centri(indeksi(i),:)); % dist je udaljenost tocke od njenog centra; d(x, c(x))
                u(i)=dist;   
                
                if(dist>l(i,j) || dist>0.5*dist_center(indeksi(i),j))
                            udaljenost=norm(A(i,:)-centri(j,:));
                            l(i,j)=udaljenost;
                            if udaljenost<dist
                                promjena=1;
                                indeksi(i)=j;
                                dist=udaljenost;
                                u(i)=udaljenost;
                            end
                        end
                    end
                end
        end
    end
    
    
    if(~promjena) break;
    end
    br_el=zeros(1,k);
    for j=1:k
        sum=zeros(1,dim);
        brojac=0;
        for i=1:n
            if(indeksi(i)==j) 
                brojac=brojac+1;
                sum=sum+A(i,:);
            end
        end
        br_el(j)=brojac;
        if brojac
            novi_centri(j,:)=sum/brojac;

        else
            br_praznih=br_praznih+1; %problem bi bio kad bi tu bio j=1;
            [maksi,maxarg]=max(br_el);
            trazi=find(indeksi==maxarg);
            ind=trazi(1);
            novi_centri(j,:)=A(ind,:);
            br_el(j)=1;
            novi_centri(maxarg,:)=(br_el(maxarg)*novi_centri(maxarg,:)-A(ind,:))/(br_el(maxarg)-1);
            br_el(maxarg)=br_el(maxarg)-1;      
            indeksi(ind)=j;
        end
    end
    
    for i=1:n
        for j=1:k
            l(i,j)=max(l(i,j)-norm(centri(j,:)-novi_centri(j,:)),0);
        end  
    end
    
    for i=1:n
        u(i)=u(i)+norm(centri(indeksi(i),:)-novi_centri(indeksi(i),:));
    end
    
    centri=novi_centri;

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