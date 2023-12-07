function entrezIDs = symbolsToEntrezIDs(Table)
%inputs:
%path to the table containing gene symbole with header symbol
   %This function takes the name of a table in which symbols are sorted
   %and returns their corresponding EntrezIDs.
   % Author: Samira Ranjbar
    T = readtable(Table);
    symbols = T.symbol;
    entrezIDs = cell(numel(symbols), 1);

    for i = 1:numel(symbols)
        symbol = symbols{i};
        entrezIDs{i} = symbolToEntrezID(symbol);
    end

    function entrezID = symbolToEntrezID(symbol)
        options = weboptions('Timeout', 20);
        url = sprintf('https://www.ncbi.nlm.nih.gov/gene/?term=%s&retmode=json', symbol);
        response = webread(url, options);
        match = regexp(response, '<span class="gene-id">ID: (\d+)</span>(.{1,200}?)\[<em>Homo sapiens<\/em> \(human\)\]', 'tokens');
        if ~isempty(match)
            entrezID = match{1}{1};
        else
            entrezID = [];
        end
    
    end
end