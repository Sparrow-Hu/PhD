function direcTest( iniTable )
%DIRECTEST Direct test on the system of equations without structural
%decomposition ; Only numerical methods are used
[conflicting,redundant]=tablerref(iniTable);
if isempty(conflicting) && isempty(redundant)
    fprintf('There is no numerical overconstraints in this system!\n');
else
    if  ~isempty(conflicting)
        fprintf('\nThere exist conflicting constraints in this system!\n');
        l=length(conflicting);
        fprintf('The number is %d,they are: \n',l);
        for i=1:l
            fprintf(' %s ', conflicting{i});
        end
        iniTable(conflicting,:)=[];
        fprintf('\nThese conflicting constraints have been removed automaticly\n');
    end
    if  ~isempty(redundant)
        fprintf('\nThere exist redundant constraints in this system\n');
        l=length(redundant);
        fprintf('The number is %d,they are: \n',l);
        for i=1:l
            fprintf(' %s ', redundant{i});
        end
        iniTable(redundant,:)=[];
        fprintf('\nThese redundant constraints have been removed automaticly\n');
    end
    
end

end

