function result = downloadAndProcessHMDBData(apiKey, hmdbFilePath, bioSamples, age)
%INPUTS:
% -exelFilePath: hmdb_ID file
% Go to https://hmdb.ca/metabolites and filter your required data, then export it.
% Make an excel file from HMDB_ID. the path to this file is excelFilePath.
% -apiKey
% To run the following code to download required concentrations you need an API key.
% Note: You need to contact the HMDB website (vasuk@ualberta.ca) to access the Human Metabolome Database (HMDB) API.
%- bioSamples: A cell array that shows the order of biosamples
      % for example:  bioSamples = {'Cerebrospinal Fluid (CSF)', 'Blood', 'Urine', 'Feces','Salvia' }
% -age:
%         possible choices:                

%     {'Newborn (0-30 days old)'     }
%     {'Infant (0-1 year old)'       }
%     {'Children (1-13 years old)'   }
%     {'Adolescent (13-18 years old)'}
%     {'Adult (>18 years old)'       }

%Author: Samira Ranjbar
%=====================================================================================================================      
   %% Step 1: Download data
   % % T is a table with HMDB ID as the first column 
T = readtable(hmdbFilePath);

% Initialize a cell array to store the response data
responseData = cell(1, 3);

options = weboptions('Timeout', 200); % Set a timeout value of 500 seconds

% Extract HMDB IDs from the table
hmdb_ids = cellstr(T{1:end, 1});

for i = 1:numel(hmdb_ids)
    hmdb_id = hmdb_ids{i};
    % Construct the URL and send the request
    url = strcat('http://35.184.189.38/api/hmdb/metabolites/concentrations/', hmdb_id, '/?api-key=', apiKey);
    
    try
        response = webread(url, options);
        
        % Process the response here
        if isfield(response, 'normal_concentration')
            % Append the data to the response cell array
            responseData{end+1, 1} = hmdb_id;
            responseData{end, 2} = response.normal_concentration;
            responseData{end, 3} = response.abnormal_concentration;
        end
    catch exception
        fprintf('Error for HMDB_ID %s: %s\n', hmdb_id, exception.message);
    end
end

% Remove the first empty row from the response cell array
responseData(1, :) = [];

% Convert the response cell array to a table
csf_mets_concentration = cell2table(responseData, 'VariableNames', {'HMDB_ID', 'normal_concentration', 'abnormal_concentration'});
%==================================================================================================

    %% Step 2: Process data
    % Perform your processing on csf_mets_concentration
    % For example, calculate the mean concentration for each metabolite
    % Initialize arrays to store the iteration index and mean values
iterationIndex = zeros(height(csf_mets_concentration), 1);
meanValues = NaN * ones(height(csf_mets_concentration), 1);
k=1;
for i = 1:height(csf_mets_concentration)
    try
        % Load your table 
        %you  need to convert structyure to a table
        T = struct2table(csf_mets_concentration.normal_concentration{i,1});
        
        % Check if the current table is empty
        if isempty(T)
            disp(['Table for iteration ', num2str(i), ' is empty. Skipping to the next iteration.']);
            continue;  % Continue to the next iteration
        end

        % Check if 'Cerebrospinal Fluid (CSF)' data exists, if not, check 'Blood', 'Urine','Feces' and then 'Salvia'
        % You can change the order of biospecimen simply by reordering
        
        % bioSamples
%         bioSamples = {'Cerebrospinal Fluid (CSF)', 'Blood', 'Urine', 'Feces','Salvia' };
        found = false;
        for j = 1:numel(bioSamples)
            
            filteredTable = T(strcmp(T.age, age) & strcmp(T.biospecimen, bioSamples{j}), :);
            if ~isempty(filteredTable)
                found = true;
                break;
            end
        end

        if found
            % Extract the 'value' column for calculation
            values = zeros(size(filteredTable, 1), 1);

            for j = 1:size(filteredTable, 1)
                % Extract numeric values from the 'value' column
                valueStr = filteredTable.value{j};

                % Check if valueStr is empty
                if isempty(valueStr)
                    % Skip this entry
                    continue;
                end
            %For different formats of concentration
                if contains(filteredTable.value{j}, '+/-')
                   numericValues = regexp(valueStr, '(\d+(\.\d+)?)', 'match','once');
                elseif contains(filteredTable.value{j}, '(')
                   numericValues = extractBefore(valueStr, '(');
                elseif contains(filteredTable.value{j}, '-')
                   numericValues = regexp(valueStr, '(\d+(\.\d+)?)', 'match');
                elseif contains(valueStr, '< ')
                numericValues = extractBetween(valueStr, '< ','uM');
                elseif contains(valueStr, '<')
                numericValues = extractBetween(valueStr, '<','uM');
                    elseif contains(filteredTable.value{j}, 'uM')
                    numericValues = extractBefore(valueStr, 'uM');
                end
               

               numericValues = str2double(numericValues);

                % Exclude NaN values (indicating missing or non-convertible values)
                numericValues = numericValues(~isnan(numericValues));

                % Store the numeric values
                values(j) = mean(numericValues, 'omitnan'); % Calculate the mean, ignoring NaN
            end

            % Calculate the mean value, ignoring NaN values
            averageValue = mean(values, 'omitnan');
        else
            % If no data found for 'Cerebrospinal Fluid (CSF)', 'Blood', 'Urine', 'Feces' or 'Salvia', set NaN as the average value
            averageValue = NaN;
        end

        % Store the iteration index and mean value
        iterationIndex(i) = i;
        meanValues(i) = averageValue;
        
    catch
        % Catch the error and display a message
        disp(['Error processing iteration ', num2str(i), '. Skipping to the next iteration.']);
        B(k)=i;
        k=k+1;
        continue;  % Continue to the next iteration
    end
end

disp(['Try again for the above metabolites...']);
Values = NaN * ones(height(csf_mets_concentration), 1);
for m = 1:length(B)
        try
        % Catch the error and display a message
            T = struct2table(csf_mets_concentration.normal_concentration{B(m),1});
            valueStr = T.value;
            if contains(valueStr, '+/-')
               numericValues = regexp(valueStr, '(\d+(\.\d+)?)', 'match','once');
            elseif contains(valueStr, '(')
               numericValues = extractBefore(valueStr, '(');
            elseif contains(valueStr, '-')
               numericValues = regexp(valueStr, '(\d+(\.\d+)?)', 'match');
            elseif contains(valueStr, '<')
                numericValues = extractBetween(valueStr, '< ','uM');
                elseif contains(valueStr, 'uM')
                numericValues = extractBefore(valueStr, 'uM');
            end
                numericValues = str2double(numericValues);
    
                % Exclude NaN values (indicating missing or non-convertible values)
                numericValues = numericValues(~isnan(numericValues));
    
                % Store the numeric values
                Values(B(m)) = mean(numericValues, 'omitnan');
        catch
             disp(['Error processing iteration ', num2str(B(m)), '. Skipping to the next iteration.']);

        end
end
%% integrate data and write the results
metMeanConcentration = NaN * ones(height(csf_mets_concentration), 1);
for i = 1: height(csf_mets_concentration) 
    if ~isnan(meanValues(i))
        metMeanConcentration(i) = meanValues(i);
    elseif isnan(meanValues(i))
        metMeanConcentration(i) = Values(i);
    end
end

 

   % Store the iteration index and mean value
        iterationIndex(i) = i;
   % Create a table with iteration index and mean values
outputTable = table(iterationIndex,csf_mets_concentration.HMDB_ID, metMeanConcentration, 'VariableNames', {'IterationIndex','HMDB_ID', 'MeanValue'});

% Write the output table to an Excel file
result = writetable(outputTable, 'output_means.xlsx');

end
