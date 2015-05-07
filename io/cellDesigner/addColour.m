function [ref] = addColour(parsed,listRxn_Color)
%
%Change colour attributes of the reaction links in a parsed CellDesigner
%model structure given a list of reaction IDs
%
%
%INPUTS
% parsed          A parsed model structure generated by 'parseCD'function
% listRxn_Color   A list of reaction IDs that need to be highilighted by
%                 changing the colour attributes of the reaciton links in
%                 the CellDesigner model. The first column stores a list
%                 of reaction IDs whose reaction links need to be
%                 highlighted, whereas the second column saves a list of
%                 Html Colours.
%
%OUTPUT
%
% ref             An updated parsed CellDesigner model structure
%
%
% Longfei Mao Oct/2014

ref=parsed;
listRxn=listRxn_Color(:,1); % coloumn 1: reaction IDs;
a=size(listRxn_Color)
if length(a)<2
    listColor=[];
    listColor=listRxn_Color(:,2); % coloumn 2: html colour codes;
    p=1;
else
    p=0;
end
for r=1:length(listRxn);    
    newRxnName=listRxn{r};    
    id=find(ismember(ref.r_info.ID,listRxn(r)))    
    if id        
        [m,n]=size(ref.r_info.ID);        
        if id>m*(n-1)  % the third column of ref.r_info.ID contains reaction ID; for example: {'re5160','re8','DESAT16_2'}
            newRxnName=ref.r_info.ID{id-2*m}; % newRxnName is defined to be the ID in the first column of ref.r_info.ID.
            column=3;
        end
    end    
    if ~column==3;                
        if ~isfield(ref,listRxn{r})
            
            newRxnName=strcat('R_',listRxn{r});
            if  isempty(strfind(newRxnName,'(e)'))
                newRxnName=strrep(newRxnName,'(e)','_e');
            end
        end
    end
    if ~isfield(ref,newRxnName)
        disp(listRxn{r});
        fprintf('error ! the listRxn{%d}',r);
        r=r+1
    else
        [rw,cw]=size(ref.(newRxnName).color)        
        for  ddr=1:rw
            if ischar(ref.(newRxnName).width{ddr,1}) % if it is char such as ('1.0');then convert it into double;
                ref.(newRxnName).width{ddr,1}=str2double(ref.(newRxnName).width{ddr,1})
            end
            %                 try
            w(1)=ref.(newRxnName).width{ddr,1}
            %             catch
            % %                 disp(w(1))
            %                  disp(newRxnName);
            %                 disp(ref.(newRxnName).width{ddr,1});
            %             end
            
            %disp(ref.(newRxnName).width(ddr,1));disp('dddd');
            %disp(w(1));
            %              if isempty(w)
            %                  w=0;
            %              end            %            
            fprintf('w value is %d\n',w(1));
            if w(1)>1&&w(1)<=10;  % flux ranges from 2 to 10 will be highligthed                
                if p==0;
                    %                     the colours (Hex triplet, e.g.,
                    %                     https://closedxml.codeplex.com/wikipage?title=Excel%20Indexed%20Colors)
                    if w(1)<=2;
                        colorStr='#FF0000FF' % flux between 1 and 2 will be highlighted in blue.
                    elseif w(1)<=5&&w(1)>2;
                        colorStr='#FFFF00FF' % flux between 2 and 5 will be highlighted in pink.
                    elseif w(1)>5;                        
                        colorStr='#FFFF0000'; % flux between 2 and 5 will be highlighted in red.
                    end                    
                elseif p==1
                    colorStr=listColor{ddr};  % ddr is the row number for each reaction;
                end                
                fprintf('p value is %d\n',p);
                colorStr=strrep(colorStr,'#','');
                fprintf('colorStr is %s\n',colorStr);                
                for  ddc=1:cw                    
                    ref.(newRxnName).color{ddr,ddc}=colorStr;
                    fprintf('set %s ''s colour to %d \n',newRxnName,ref.(newRxnName).color{ddr,ddc});
                end
            end
        end        
    end    
end

