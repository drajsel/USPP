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
% u centrima je postavljena pocetna inicijalizacija
br_iter=1; % u ostalim kodovima jedno pridruzivanje klasterima je prije petlje
tic
br_praznih=0;

promjena=1;
br_iter=0; 
indeksi=zeros(1,n);
while promjena
    promjena=0;
    br_iter=br_iter+1;
    
    for i=1:n
        dist=zeros(1,k); % udaljenosti od tocke x do centra j
        for j=1:k
            dist(j)=norm(A(i,:)-centri(j,:));
        end
        [mini, ind]=min(dist); % ind je vektor za koji se poprima minimalna vrijednost
        if br_iter==1
            indeksi(i)=ind;
            promjena=1;
        elseif mini<dist(indeksi(i))
                indeksi(i)=ind;
                promjena=1;
        end
             
    end
    % sad radimo nove centroide koji su aritmetièka sredina svih vektora koji mu pripadaju
	stari=centri;
	centri=zeros(k,dim);
	br_el=zeros(1,k);
for j=1:k
	trazi=find(indeksi==j);
	br_el(j)=size(trazi,2);
	if ~br_el(j)
		br_praznih=br_praznih+1;
		[maksi,maxarg]=max(br_el);	
		trazi=find(indeksi==maxarg);
            ind=trazi(1);
            centri(j,:)=A(ind,:);
            br_el(j)=1;
            centri(maxarg,:)=(br_el(maxarg)*centri(maxarg,:)-A(ind,:))/(br_el(maxarg)-1);
            br_el(maxarg)=br_el(maxarg)-1;      
            indeksi(ind)=j;
    else
    for i=1:size(trazi,2)
        centri(j,:)=centri(j,:)+A(trazi(i),:);
    end
    centri(j,:)=centri(j,:)/size(trazi,2);
    end
end
end
toc
display(br_iter);

x=input('Zelite li plot? \n 0: ne\n 1: da\n');
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