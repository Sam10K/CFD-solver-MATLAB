function[BC] = bc_read()
BC = struct;

for file = ["U.bc","V.bc","P.bc"]
    
    fid = fopen(strcat('BC/',file),'r');
    field = split(file,'.bc');
    field = field{1};
    data = fgetl(fid);

    while ~feof(fid)

        dummy = sscanf(data,'ZONE-%s');
        if ~isempty(dummy)
            zone = dummy;
        end

        dummy = sscanf(data,' type %s;');
        if ~isempty(dummy)
            type = split(dummy,';');
            type = type{1};
            BC.(zone).(field).type = type;
        end

        dummy = sscanf(data,' value %s;');
        if ~isempty(dummy)
            value = split(dummy,';');
            value = value{1};
            BC.(zone).(field).value = str2double(value);
        end
        
        dummy = sscanf(data,' list_type %s;');
        if ~isempty(dummy)
            list_type = split(dummy,';');
            list_type = list_type{1};
            if strcmp(list_type,'uniform')
                curr = fgetl(fid);
                list_values = sscanf(curr,' list_values %g;');
            elseif strcmp(list_type,'non-uniform')
                curr = fgetl(fid);
                nval = sscanf(curr,' list_values %d(');
                list_values = zeros(nval,1);
                i=1;
                curr = fgetl(fid);
                while ~strcmp(curr,');')
                    list_values(i)=sscanf(curr,'%g');
                    curr = fgetl(fid);
                    i=i+1;
                end
            end  
        end
            
        data = fgetl(fid);
    end
    
    fclose(fid);
    BC.internal.(field).value = list_values;    
    
end

end