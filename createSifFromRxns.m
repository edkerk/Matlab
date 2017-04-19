function createSifFromRxns(model,rxns,out,compartments,use_rev,filename,isCobra,maxDegree)

% function createSifFromRxns
%
% This function creates a Cytoscape-readable .sif file for
% metabolite-gene/reaction networks. Only a GEM must be specified. If run
% on default, the resulting file is a non-compartmentalized metabolite-gene 
% non-directional network of all genes found in the model.
% The output is saved in the current directory with the name
% MetsToGenes.sif'.
%
% INPUT
% model:    RAVEN-format GEM
% rxns:     (optional) A cell array of strings containing the list of
%           reactions for which the network should be reconstructed. Default is all
%           reactions in the model. The cell array of strings must match the .rxns
%           field in the model to be included in the file.
% out:      (optional) A string, either 'rm','gm' (default), 'rg', 'mm','rr' or 'pg', whether a
%           metabolite-reaction network, a metabolite-gene network, a
%           gene-reaction, a metabolite-metabolite, a rxn-rxn, or a gene-pathway network should be created.
% compartments: (optional) a logical, indicating whether metabolites in the
%               file should be considered unique if belonging to different compartments.
%               Default is false.
% use_rev:  (optional) A logical, indicating whether reversible reactions are 
%           written in the file in both directions. Default is false.
% filename: (optional) A string for the output filename. Default is
%           'MetsToGenes.sif'
% isCobra: (optional) A logical, if the model is in COBRA-format (def: false)
% maxDegree: (optional) Only if out is "rr". The maximum degree for a
%           metabolite to declare the rxn involving the metabolite as
%           connected (def:Inf).
%
% USAGE
%
% createSifFromRxns(model,rxns,out,compartments,use_rev,filename,isCobra,maxDegree)
%
% Francesco Gatto - 2014-11-11

if nargin < 8 || isempty(maxDegree)
    maxDegree = Inf;
end
if nargin < 7 || isempty(isCobra)
    isCobra = false;
end
if nargin < 6 || isempty(filename)
    fid = fopen('MetsToGenes.sif','wt');
    fprintf('\nFile will be saved in "MetsToGenes.sif".\n')
else
    fid = fopen(filename,'wt');
end
if nargin < 5 || isempty(use_rev)
    use_rev = false;
end
if nargin < 4 || isempty(compartments)
    compartments = false;
end
if nargin < 3 || isempty(out)
    out = 'gm';
end
if nargin < 2 || isempty(rxns)
    rxns = model.rxns;
end

if ~strcmp(out,{'rm';'gm';'rg';'pg';'mm';'rr'})
    fprintf('\nError: invalid out argument. It must be either "rm", "gm", "rg", "rr", or "pg".\n')
    return
end

genes     = model.genes;
if ~isCobra
    mets      = model.metNames;
else
    mets = model.mets;
end
rev       = model.rev;
nMets     = length(mets);

if ~isCobra
    if compartments 
        comps       = model.comps(model.metComps); %model.compNames(str2double(model.metComps));
        metNamesC   = cellfun(@(a,b,c,d) [a,b,c,d],mets,repmat({'['},nMets,1),comps,repmat({']'},nMets,1),'uni',false);
        metNames    = metNamesC;
    else
        metNames    = mets;
    end
else
    if compartments
        metNames = mets;
    else
        fprintf('\Warning: not compartmentalized mets are not supported in COBRA models.\n')
    end
end

nRxns = length(rxns);
unmappedCounter = 0;
if ~strcmp('rr',out)
    for i = 1:nRxns
        currentRxn = rxns(i);
        [mapped,indRinS] = ismember(currentRxn,model.rxns);
        if ~mapped
            unmappedCounter = unmappedCounter + 1;
            unmappedRxns(unmappedCounter) = currentRxn;
            continue
        else
            indProductsinRxn    = model.S(:,indRinS)>0;
            indSubstratesinRxn  = model.S(:,indRinS)<0;
            productsInRxn     	= metNames(indProductsinRxn);
            substratesInRxn     = metNames(indSubstratesinRxn);
            nSubs = length(substratesInRxn);
            nProd = length(productsInRxn);
            switch out
                case 'rm'
                    switch use_rev                    
                        case 0
                            for j = 1:nSubs
                                fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',currentRxn{:});
                            end
                            for k = 1:nProd
                                fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',productsInRxn{k});
                            end
                        case 1
                            if rev(indRinS)
                                for j = 1:nSubs
                                    fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',currentRxn{:});
                                    fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',substratesInRxn{j});
                                end
                                for k = 1:nProd
                                    fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',productsInRxn{k});
                                    fprintf(fid,'%s\t%s\t%s\n',productsInRxn{k},'mr',currentRxn{:});
                                end
                            else
                                for j = 1:nSubs
                                    fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mr',currentRxn{:});
                                end
                                for k = 1:nProd
                                    fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'mr',productsInRxn{k});
                                end
                            end
                    end
                case 'gm'
                    indGenesinRxn = model.rxnGeneMat(indRinS,:)~=0;
                    genesInRxn    = genes(indGenesinRxn);
                    nGenes        = length(genesInRxn);
                        switch use_rev                    
                            case 0
                                for l = 1:nGenes
                                    for j = 1:nSubs
                                        fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'gm',genesInRxn{l});
                                    end
                                    for k = 1:nProd
                                        fprintf(fid,'%s\t%s\t%s\n',productsInRxn{k},'gm',genesInRxn{l});
                                    end
                                end
                            case 1
                                if rev(indRinS)
                                    for l = 1:nGenes
                                        for j = 1:nSubs
                                            fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'gm',genesInRxn{l});
                                            fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'gm',substratesInRxn{j});
                                        end
                                        for k = 1:nProd
                                            fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'gm',productsInRxn{k});
                                            fprintf(fid,'%s\t%s\t%s\n',productsInRxn{k},'gm',genesInRxn{l});
                                        end
                                    end
                                else
                                    for l = 1:nGenes
                                        for j = 1:nSubs
                                            fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'gm',genesInRxn{l});
                                        end
                                        for k = 1:nProd
                                            fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'gm',productsInRxn{k});
                                        end
                                    end
                                end
                        end
                case 'pg'
                    indGenesinRxn = model.rxnGeneMat(indRinS,:)~=0;
                    genesInRxn    = genes(indGenesinRxn);
                    nGenes        = length(genesInRxn);
                    currentSubs   = model.subSystems(indRinS);
                        switch use_rev                    
                            case 0
                                for l = 1:nGenes
                                    for j = 1:nSubs
                                        fprintf(fid,'%s\t%s\t%s\n',currentSubs{:},'pg',genesInRxn{l});
                                    end
                                    for k = 1:nProd
                                        fprintf(fid,'%s\t%s\t%s\n',currentSubs{:},'pg',genesInRxn{l});
                                    end
                                end
                            case 1
                                if rev(indRinS)
                                    for l = 1:nGenes
                                        for j = 1:nSubs
                                            fprintf(fid,'%s\t%s\t%s\n',currentSubs{:},'pg',genesInRxn{l});
                                            fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'pg',currentSubs{:});
                                        end
                                        for k = 1:nProd
                                            fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'pg',currentSubs{:});
                                            fprintf(fid,'%s\t%s\t%s\n',currentSubs{:},'pg',genesInRxn{l});
                                        end
                                    end
                                else
                                    for l = 1:nGenes
                                        for j = 1:nSubs
                                            fprintf(fid,'%s\t%s\t%s\n',currentSubs{:},'pg',genesInRxn{l});
                                        end
                                        for k = 1:nProd
                                            fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'pg',currentSubs{:});
                                        end
                                    end
                                end
                        end
                case 'rg'
                indGenesinRxn = model.rxnGeneMat(indRinS,:)~=0;
                genesInRxn    = genes(indGenesinRxn);
                nGenes        = length(genesInRxn);
                    switch use_rev                    
                        case 0
                            for l = 1:nGenes
                                fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'rg',genesInRxn{l});
                            end
                        case 1
                            if rev(indRinS)
                                for l = 1:nGenes
                                    fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'rg',genesInRxn{l});
                                    fprintf(fid,'%s\t%s\t%s\n',genesInRxn{l},'rg',currentRxn{:});
                                end
                            else
                                for l = 1:nGenes
                                    fprintf(fid,'%s\t%s\t%s\n',currentRxn{:},'rg',genesInRxn{l});
                                end
                            end
                    end
                case 'mm'
                switch use_rev                    
                    case 0
                        for j = 1:nSubs
                            for k = 1:nProd
                                fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mm',productsInRxn{k});
                            end
                        end
                    case 1
                        if rev(indRinS)
                            for j = 1:nSubs
                                for k = 1:nProd
                                    fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mm',productsInRxn{k});
                                    fprintf(fid,'%s\t%s\t%s\n',productsInRxn{k},'mm',substratesInRxn{j});
                                end
                            end
                        else
                            for j = 1:nSubs
                                for k = 1:nProd
                                    fprintf(fid,'%s\t%s\t%s\n',substratesInRxn{j},'mm',productsInRxn{k});
                                end
                            end
                        end
                end

            end
        end
    end
else
    for i = 1:length(mets)
        interactRxnInd = find(model.S(i,:));
        if numel(interactRxnInd)>1 && numel(interactRxnInd)<maxDegree%Skip dead-end rxns or high degree mets
            interactRxnPairs = nchoosek(interactRxnInd,2);
            for pair = 1:size(interactRxnPairs,1)
                fprintf(fid,'%s\t%s\n',model.rxns{interactRxnPairs(pair,1)},model.rxns{interactRxnPairs(pair,2)});
            end
        end
    end    
end
fclose(fid);

if unmappedCounter > 0
    fprintf('\nThe following reactions could not be matched in the input model: ')
    disp(unmappedRxns)
    fprintf('\n')
else
    fprintf('\nAll reactions were matched to input model.\n')
end