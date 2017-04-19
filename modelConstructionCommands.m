% As you will see, we can add multiple reactions at the same time, and the same goes for genes, metabolites etc.
% However, it will become difficult to track all the changes if you make many changes all at once. It's therefore
% better to perhaps add a comment in your script, giving a short explanation why a certain change is made,
% followed by a series of comments specific about that change. For instance: these genes are substituted as
% manually found reciprocal blast hits; these reactions are added to include this particular pathway, etc.

% How to change the gene association for a reaction: Make sure you change the model name,
% reaction ID and gene associations. Use 'and' and 'or' instead of ':' and ';'. Use brackets ()
% to specify priority of certain relationships over others (as done below). With this command
% you overwrite the previous gene assocation for this reaction only. If genese should be deleted
% from the old gene association, just don't mention it in the new gene association. If genes occur
% in multiple reactions, make sure you modify each reaction individually.

% The changeGeneAssoc command is not part of RAVEN, make sure you have it in a Matlab path (or your
% current working directory).

model=changeGeneAssoc(model,'y000001','(YALI0C06446g and YALI0D09273g) or (YALI0E03212g and YALI0D09273g)',true);

% changeGeneAssoc will also automatically add new genes that had not been used in the model.
% However, you will have to manually remove old genes that you no longer want to include. With
% removeGenes you can remove multiple genes at the same time:

model=removeGenes(model,{'YALI0B14619g','OLD_sce_YAL150W'},false,false);

% If you want to add new reactions, you need to first specify all the information you have for the 
% reactions (here we're adding two reactions:):

rxnToAdd.rxns={'y200066','y200067'}; % These are new reactions, so you need to specify reaction IDs yourself.
rxnToAdd.rxnNames={'erythrose kinase',...
	'erythrose reductase'}; % Each reaction also needs a name, that is more descriptive than the reaction ID
	% ... can be used after a comma if you want to split a command over multiple lines, to make
	% it easier to read
rxnToAdd.eccodes={'','1.1.1.21'}; % EC numbers are not compulsory, but are good to add if you know them.
rxnToAdd.subSystems={'Pentose and glucuronate interconversions',...
	'Pentose and glucuronate interconversions'}; % The subsystem describes what pathway the reactions belong to.
	% Can be an exisiting subSystem, but you can also choose to introduce a new subsystem specifically for
	% your new reactions.
rxnToAdd.grRules={'(YALI0C06446g and YALI0D09273g)','YALI0B16192g'}; % grRules as described above
rxnToAdd.equations={'D-erythrose 4-phosphate[cy] + ADP[cy] <=> erythrose[cy] + ATP[cy]',...
	'erythrose[cy] + NADH[cy] + H+[cy] => erythritol[cy] + NAD[cy]'}; % The new reactions are written
	% out in full. Pay attention to adding the compartments, using the correct metabolite names
	% and the directionality of the reaction (indicated by the arrow)
genesToAdd.genes={'YALI0B16192g'}; % If there are genes in the grRules that had not been used before you
% should specify and add them separately.
model=addGenes(model,genesToAdd);
model=addRxns(model,rxnToAdd,3,{},true); % This command will add the new reactions.

% If you're adding new reactions, you probably also add new metabolites. Using the approach above, new
% metabolites will be automatically added, but there will not be much information included in the model
% for these metabolites. If you do have more information, you can add it with the following commands, but
% make sure you run this BEFORE the addRxns command described above.

metsToAdd.metNames={'3-methylbutanoyl-CoA','2-methylpropanoyl-CoA'}; % Metabolite names
metsToAdd.compartments={'mi','mi'}; % Abbreviation of which compartment they should reside
metsToAdd.metFormulas={'C26H44N7O17P3S','C25H42N7O17P3S'}; % Chemical formula of each metabolite.
metsToAdd.inchis={'1S/C26H44N7O17P3S/c1-14(2)9-17(35)54-8-7-28-16(34)5-6-29-24(38)21(37)26(3,4)11-47-53(44,45)50-52(42,43)46-10-15-20(49-51(39,40)41)19(36)25(48-15)33-13-32-18-22(27)30-12-31-23(18)33/h12-15,19-21,25,36-37H,5-11H2,1-4H3,(H,28,34)(H,29,38)(H,42,43)(H,44,45)(H2,27,30,31)(H2,39,40,41)/t15-,19-,20-,21+,25-/m1/s1',...
'1S/C25H42N7O17P3S/c1-13(2)24(37)53-8-7-27-15(33)5-6-28-22(36)19(35)25(3,4)10-46-52(43,44)49-51(41,42)45-9-14-18(48-50(38,39)40)17(34)23(47-14)32-12-31-16-20(26)29-11-30-21(16)32/h11-14,17-19,23,34-35H,5-10H2,1-4H3,(H,27,33)(H,28,36)(H,41,42)(H,43,44)(H2,26,29,30)(H2,38,39,40)/t14-,17-,18-,19+,23-/m1/s1'}; % InChis are accurate descriptions of the chemical
% structure of the metabolite, much more accurate than the formula. InChis can for instance be found
% on ChemSpider (http://www.chemspider.com/)
metsToAdd.metMiriams={'obo.chebi:CHEBI:15487','obo.chebi:CHEBI:15479'}; % Miriams are links to other
% databases. Here there is a mention of the metabolite ID in the ChEBI database
% (https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:15487)
model=addMets(model,metsToAdd); % Adds the metabolites to the model. Should be done BEFORE addRxns, if
% you want this more detailed information about the new metabolites included in the model (otherwise addRxns
% will have already added the metabolites with no additional information).

% If you are certain that a particular reaction should be removed, you can use the following command:
model=removeRxns(model,{'1542','20'},true,true,true); % Where the numbers here correspond to reaction IDs.

% If you want to change a reaction, because something was wrong:
model=changeRxns(model,'y001076',{'ureidoglycolic acid[cy] => urea[cy] + glyoxylate[cy]'},3); % This changes
% the previous reaction equation with the new reaction equation.



