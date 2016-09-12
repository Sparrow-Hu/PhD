function [ Conflicting,Redundant ] = qrRC( table,compDAG,tol )
%QRRC Identify the redundant  and conflicting constraints based on qr
%factorization;theoritically it is the most stable one comparing to
%LU/Gauss Jordan Elimination
Conflicting={};
Redundant={};
A=table2array(table(:,1:end-1));  
%% transpose the A for column manipulation in qr factorization
A=A';
[m,n]=size(A);
[Q,R,E] = qr(A,0);% R is the new coordinate matrix of new orthogonal basis: that is, A is projected onto Q column space spanned by its columns which are orthogonal
table=table(E,:);%row permutation so that singular values are in decreasing order
B=table2array(table(:,end)); 
B=B';
A=A(:,E);
r=rank(A);
[mm,nn]=size(R);
if r< n
    temp1=1;
    temp2=1;
    relationMatrix=inv(R(1:r,1:r))*R(1:r,r+1:end);% A(:,r+1:end)=A(:,1:r)*relationMatrix
    residue_B=B(r+1:end)-B(1:r)*relationMatrix;% find the difference between b values
    residue_B=abs(residue_B);
    for i=r+1:nn
        group={};
            t1=1;
            for j=1:r 
                if abs(relationMatrix(j,i-r))>tol
                    group(t1)=table.Properties.RowNames(j);% find the basis equations that form linear conbination of the previous conflicting equation
                    t1=t1+1;
                end
            end
            source_group=findSource(compDAG,group);
        if residue_B(i-r)> tol 
            Conflicting(temp1) = table.Properties.RowNames(i);%equations that are not satisfied
            source_C=findSource(compDAG,Conflicting(temp1));
            [m,n]=size(source_group);
            if m ~= 1
                source_group = source_group';
            end
            [M,N]=size(source_C);
            if M ~= 1
                source_C = source_C';
            end
            Group=unique([group,source_group,source_C]);
            if ~isempty(Group)
                disp([Conflicting(temp1),'is conflicting with group:',Group]);
            end
            temp1=temp1+1;
        else
            Redundant(temp2) = table.Properties.RowNames(i);% equations that are satisfied but are redundant
            source_R=findSource(compDAG,Redundant(temp2));
            [m,n]=size(source_group);
            if m ~= 1
                source_group = source_group';
            end
            [M,N]=size(source_R);
            if M ~= 1
                source_R = source_R';
            end
            Group=unique([group,source_group,source_R]);
            if ~isempty(Group)
                disp([Redundant(temp2),'is redundant with group:',Group]);
            end
            temp2=temp2+1;
        end
    end
end
end

