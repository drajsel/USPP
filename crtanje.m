

for i=1:n
    plot(A(i,1), A(i,2), '.b');
    if(i==1) hold on;
    end
end

for i=1:k
    plot(centri(i,1), centri(i,2), '+r');
end


hold off