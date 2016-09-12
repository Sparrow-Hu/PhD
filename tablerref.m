function [Conflicting,Redundant] = tablerref( T,tol,type )
%TABLERREF Find the basis of equation space and output the overconstraints
%   this function is  based on the frref function in the package TEST
Conflicting={};
Redundant={};
A=table2array(T);
[m,n]=size(A);

P=linspace(1,m,m);
% A=A';
% T=array2table(A,'VariableNames',T.Properties.RowNames,'RowNames',T.Properties.VariableNames);
%   R = FRREF(A) produces the reduced row echelon form of A.
%   [R,jb] = FRREF(A,TOL) uses the given tolerance in the rank tests.
%   [R,jb] = FRREF(A,TOL,TYPE) forces frref calculation using the algorithm
%   for full (TYPE='f') or sparse (TYPE='s') matrices.
%   
% 
%   Description: 
%   For full matrices, the algorithm is based on the vectorization of MATLAB's
%   RREF function. A typical speed-up range is about 2-4 times of 
%   the MATLAB's RREF function. However, the actual speed-up depends on the 
%   size of A. The speed-up is quite considerable if the number of columns in
%   A is considerably larger than the number of its rows or when A is not dense.
%
%   For sparse matrices, the algorithm ignores the tol value and uses sparse
%   QR to compute the rref form, improving the speed by a few orders of 
%   magnitude.
%
%   Authors: Armin Ataei-Esfahani (2008)
%            Ashish Myles (2012)
%
%   Revisions:
%   25-Sep-2008   Created Function
%   21-Nov-2012   Added faster algorithm for sparse matrices

[m,n] = size(A);

switch nargin
  case 1,
    % Compute the default tolerance if none was provided.
%     tol = max(m,n-1)*eps(class(A))*norm(A(:,1:end-1),'inf');
     tol = max(m,n-1)*eps(norm(A(:,1:end-1),'inf'));
    if issparse(A)
      type = 's';
    else
      type = 'f';
    end
  case 2,
  if issparse(A)
    type = 's';
  else
    type = 'f';
  end
  case 3,
    if ~ischar(type)
      error('Unknown matrix TYPE! Use ''f'' for full and ''s'' for sparse matrices.')
    end
    type = lower(type);
    if ~strcmp(type,'f') && ~strcmp(type,'s')
      error('Unknown matrix TYPE! Use ''f'' for full and ''s'' for sparse matrices.')
    end
end

r=rank(A(:,1:end-1));%add by Hao HU

do_full = ~issparse(A) || strcmp(type,'f');

if do_full
    % Loop over the entire matrix.
    i = 1;
    j = 1;
    jb = [];
    % t1 = clock;
    while (i <= m) && (j < n)
       % Find value and index of largest element in the remainder of column j.
       [p,k] = max(abs(A(i:m,j))); k = k+i-1;
       if (p <= tol)
          % The column is negligible, zero it out.
          A(i:m,j) = 0; %(faster for sparse) %zeros(m-i+1,1);
          j = j + 1;
       else
          % Remember column index
          jb = [jb j];
          % Swap i-th and k-th rows.
          A([i k],j:n) = A([k i],j:n);
          P([i k])=P([k i]);
          % Divide the pivot row by the pivot element.
          Ai = A(i,j:n)/A(i,j);    
          % Subtract multiples of the pivot row from all the other rows.
          A(:,j:n) = A(:,j:n) - A(:,j)*Ai;
          A(i,j:n) = Ai;
          i = i + 1;
          j = j + 1;
       end
    end
else
    % Non-pivoted Q-less QR decomposition computed by Matlab actually
    % produces the right structure (similar to rref) to identify independent
    % columns.
    R = qr(A);

    % i_dep = pivot columns = dependent variables
    %       = left-most non-zero column (if any) in each row
    % indep_rows (binary vector) = non-zero rows of R
    [indep_rows, i_dep] = max(R ~= 0, [], 2);
    indep_rows = full(indep_rows); % probably more efficient
    i_dep = i_dep(indep_rows);
    i_indep = setdiff(1:n, i_dep);

    % solve R(indep_rows, i_dep) x = R(indep_rows, i_indep)
    %   to eliminate all the i_dep columns
    %   (i.e. we want F(indep_rows, i_dep) = Identity)
    F = sparse([],[],[], m, n);
    F(indep_rows, i_indep) = R(indep_rows, i_dep) \ R(indep_rows, i_indep);
    F(indep_rows, i_dep) = speye(length(i_dep));

    % result
    A = F;
    jb = i_dep;
end
% AA=array2table(A,'VariableNames',T.Properties.VariableNames,'RowNames',T.Properties.RowNames);
% r=rank(A);
% [m,n]=size(AA);
% rr=1;
% temp=0;
% mak=find(A(:,1));
% mak=mak(end);
% for i=2:n
%     K=find(A(:,i));
%     if mak<K(end)
%         mak=K(end);
%     else
%         temp=temp+1;
%         Overconstraints(temp)=AA.Properties.VariableNames(i);
%     end
% end
T=T(P,:);%after the rref process, perturbate the rows in order to  match b which is A(:,end)
%% it is not a good idea to use tolerance as a criteria for evaluating redundant and conflicting constraints in 
temp1=0;
temp2=0;
for i=r+1:m
    if abs(A(i,end)) >=tol
        temp1=temp1+1;
        Conflicting(temp1) = T.Properties.RowNames(i);
    else
        temp2=temp2+1;
        Redundant(temp2) = T.Properties.RowNames(i);
    end
end
        

    




end

