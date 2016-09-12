fprintf(' Adding equation names and system of coefficients to adjacent table\n');
global iniTable ceqPara ceqVars 
if ~isempty(ceqPara) && ~isempty(ceqVars)
    ceqPara=vpa(ceqPara);ceqVars=vpa(ceqVars);
    % assume(ceqPara,'real');assume(ceqVars,'real');
    ceqPara=array2table(ceqPara.');
    ceqVars=array2table(ceqVars.');
end
[m,n]=size(iniTable);%compute the number of rows and columns of initial table 
rownames =[];
Rownames = {};
%% generate names for the rows(equations) of initial table
for i=1:m
    rownames = [rownames,sym(['e',num2str(i)])];
end
for i=1:m
    Rownames{i} = char(rownames(i));
end
iniTable.Properties.RowNames = Rownames;
if ~isempty(ceqPara) && ~isempty(ceqVars)
    ceqPara.Properties.RowNames = Rownames;
    ceqVars.Properties.RowNames = Rownames;
end
%% Check if there exist rows whose coefficients are all zeros
iniMatrix = table2array(iniTable);
[m,n]=size(iniMatrix);
index0s=[];
index0sC=[];
for i=1:m
    if all(iniMatrix(i,1:end-1)==0)
        index0s=[index0s,i];
        if iniMatrix(i,end)~=0
            index0sC=[index0sC,i];
        end
    end
end
if ~isempty(index0s)
    disp([iniTable.Properties.RowNames{index0s}, ' have all 0s coefficients, thus have been removed!\n']);
    disp([iniTable.Properties.RowNames{index0sC}, ' have all 0s coefficients with non-0s on the right side of equations\n']);
    iniTable(index0s,:) = [];
end
    
%% get the coefficient table T by removing the "b" column of iniTable 
T= iniTable;
[m,n]=size(iniTable);%size of new 'iniTable'
T.b=[];
%% Turn T to adjacent table where only 0 and 1 exists
for i=1:m
    for j=1:n-1
        if T{i,j}~= 0
            T{i,j}=1;
        end
    end
end
strAnalysis(T);
% %% Try LU on the entire system
% fprintf('\n Apply LU directly on the system!\n');
% [conflicting,redundant]=findRC(TOriginal,1e-4);
% %%Output conflicting/redundant constraints given by LU
% if  ~isempty(conflicting)
%     fprintf('There exist conflicting constraints:\n');
%     l=length(conflicting);
%     fprintf('The number is %d,they are: \n',l);
%     for i=1:l
%         fprintf(' %s ', conflicting{i});
%     end
% end
% if  ~isempty(redundant)
%     fprintf('There exist redundant constraints:\n');
%     l=length(redundant);
%     fprintf('The number is %d,they are: \n',l);
%     for i=1:l
%         fprintf(' %s ', redundant{i});
%     end
% end
% %% Try Gauss-Jordan Elimination on the entire system
% fprintf('\n Apply GE directly on the system!\n');
% [conflicting,redundant]=tablerref(TOriginal,1e-4);
% %%Output conflicting/redundant constraints given by LU
% if  ~isempty(conflicting)
%     fprintf('There exist conflicting constraints:\n');
%     l=length(conflicting);
%     fprintf('The number is %d,they are: \n',l);
%     for i=1:l
%         fprintf(' %s ', conflicting{i});
%     end
% end
% if  ~isempty(redundant)
%     fprintf('There exist redundant constraints:\n');
%     l=length(redundant);
%     fprintf('The number is %d,they are: \n',l);
%     for i=1:l
%         fprintf(' %s ', redundant{i});
%     end
% end
% 
% 
