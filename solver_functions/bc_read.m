%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to read the boundary condition files and store it as a structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


        data = fgetl(fid);
    end
    fclose(fid);
    
end

end