% hybrid.m
% Just input a matrix A that is n x n+1 and the number of equations n.
function [x]=hybrid(n, A)
NROW = (1:n)';	        % pivoting vector NROW
for i = 1:(n-1)
  s = max(abs(A(:,1:n)')); % max of each row row
  LHS = abs(A(NROW(i),i)/s(NROW(i)));
  kp = i;
  for j = (i+1):n
    RHS = abs(A(NROW(j),i)/s(NROW(j)));
    if RHS > LHS,  LHS = RHS;  kp = j;  end
  end
  l = NROW(kp);  NROW(kp) = NROW(i);  NROW(i) = l; 
  for j = (i+1):n
      m(NROW(j),i) = A(NROW(j),i)/A(NROW(i),i);
      A(NROW(j),:) = A(NROW(j),:)-m(NROW(j),i)*A(NROW(i),:);
  end
end

x(n) = A(NROW(n),n+1)/A(NROW(n),n);
for i=n-1:-1:1
    x(i) = A(NROW(i),n+1);
    for j=(i+1):n
        x(i) = x(i) - A(NROW(i),j)*x(j);
    end
    x(i) = x(i) / A(NROW(i),i);
end
x % solutions
end

% test.m %
A = [3.333, 15920, 10.333, 7953;
	2.222, 16.170, 9.612, 0.965;
	-1.5611, 5.1792, -1.6855, 2.714];
hybrid(3,A);



    
    
    
