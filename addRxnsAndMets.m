function newModel=addsRxnsGenesMets(model,sourcemodel,rxns,rxnNote,confidence)
% addRxnsGenesMets
%   Adds reactions to a model, including new metabolites. Does not take
%   genes into account.
%
%   model           draft model where reactions should be added
%   sourcemodel     model where reactions and metabolites are sourced from
%   rxns            cell array with reaction IDs (from source model)
%   rxnNote         string explaining why reactions were added to model,
%                   is included as newModel.rxnNotes (opt, default
%                   'Added via addRxnsAndMets()')
%   confidence      string, specifying confidence score for all reactions.
%                   Following doi:10.1038/nprot.2009.203 (opt, default 0)
%                   4:  biochemical data: direct evidence from enzymes
%                       assays
%                   3:  genetic data: knockout/-in or overexpression
%                       analysis
%                   2:  physiological data: indirect evidence, e.g.
%                       secretion products or defined medium requirement
%                       sequence data: genome annotation
%                   1:  modeling data: required for functional model,
%                       hypothetical reaction
%                   0:  no evidence
%
%   newModel         an updated model structure
%
%   Usage: newModel=addRxns(model,rxnsToAdd,eqnType,compartment,allowNewMets)
%
%   Eduard Kerkhoven, 2016-12-20

if nargin<5
    confidence=0;
end
if nargin<4
    rxnNote='Added via addRxnsAndMets()';
end


%% Obtain indexes of reactions in source model
notNewRxn=rxns(ismember(rxns,model.rxns));
rxns=rxns(~ismember(rxns,model.rxns));
if isempty(rxns)
    throw(MException('','All reactions are already in the model.'));
end

if ~isempty(notNewRxn)
    disp('The following reactions were already present in the model and will not be added:')
    disp(strjoin(notNewRxn,'\n'))
end


rxnIdx=find(ismember(sourcemodel.rxns,rxns)); % Get rxnIDs

%% Add new metabolites
metIdx=find(any(sourcemodel.S(:,rxnIdx),2)); % Get metabolite IDs
% Many of the metabolites in are already in the draft model, so only add the new metabolites
mets=sourcemodel.mets(metIdx);
notNewMet=mets(ismember(mets,model.mets));
if ~isempty(notNewMet)
    disp('The following metabolites were already present in the model and will not be added:')
    disp(strjoin(notNewMet,'\n'))
end

metIdx=metIdx(~ismember(sourcemodel.mets(metIdx),model.mets));

if ~isempty(metIdx)
    metsToAdd.mets=sourcemodel.mets(metIdx);
    metsToAdd.metNames=sourcemodel.metNames(metIdx);
    metsToAdd.metFormulas=sourcemodel.metFormulas(metIdx);

    metsToAdd.compartments=strtrim(cellstr(num2str(sourcemodel.metComps(metIdx)))); % Convert from compartment string to comparment number
    [~,idx]=ismember(metsToAdd.compartments,strsplit(num2str(1:length(sourcemodel.comps)))); % Match compartment number to compartment abbreviation
    metsToAdd.compartments=sourcemodel.comps(idx); % Fill in compartment abbreviations

    model=addMets(model,metsToAdd);
end
disp('Number of metabolites added to the model:')
disp(numel(metIdx))

%% Add new genes
rxnToAdd.grRules=sourcemodel.grRules(rxnIdx);
geneList=regexp(rxnToAdd.grRules,'[)(]*|( and )*|( or )*','split');
genesToAdd.genes=setdiff(unique(horzcat(geneList{:})),sourcemodel.genes);
model=addGenes(model,genesToAdd);

disp('Number of genes added to the model:')
disp(numel(genesToAdd.genes))
%% Add new reactions
rxnToAdd.equations=constructEquations(sourcemodel,rxnIdx);
rxnToAdd.rxnNames=sourcemodel.rxnNames(rxnIdx);
rxnToAdd.rxns=sourcemodel.rxns(rxnIdx);
rxnToAdd.grRules=sourcemodel.grRules(rxnIdx);
rxnToAdd.lb=sourcemodel.lb(rxnIdx);
rxnToAdd.ub=sourcemodel.ub(rxnIdx);
rxnToAdd.rxnNotes=cell(1,numel(rxnToAdd.rxns));
rxnToAdd.rxnNotes(:)={rxnNote};
rxnToAdd.confidenceScores=cell(1,numel(rxnToAdd.rxns));
rxnToAdd.confidenceScores(:)={confidence};
newModel=addRxns(model,rxnToAdd,3,'',false);

disp('Number of reactions added to the model:')
disp(numel(rxnIdx))

