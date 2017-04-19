function [ rxngene ] = getRxnGeneList(model)

n=0;
rxngene=cell(2,1);
for i=1:length(model.rxns)
    for j=1:length(model.genes)
        if model.rxnGeneMat(i,j)==1
            n=n+1;
            rxngene(n,1)=model.rxns(i);
            rxngene(n,2)=model.genes(j);
        end
    end
end
clear i j n
end
