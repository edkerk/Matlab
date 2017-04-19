cd('D:\Users\eduardk.NET\Dropbox\Postdoc\Annotate YLI');
[~,~,GOs]=xlsread('blast2go_table_20130501_0947.csv');
l=0
%output=zeros(length(GOs),2);
for i=1:length(GOs) % Go through all genes
    for j=2:size(GOs,2) % Go through all GO terms per gene
        if ~isnan(GOs{i,j})
            if ~strcmp(GOs(i,j),'-')
                l=l+1 % Row number for output file
                output(l,1)=GOs(i,1);
                output(l,2)=GOs(i,j);
            end
        end
    end
end
xlswrite('temp.xlsx',output)    
        
        
