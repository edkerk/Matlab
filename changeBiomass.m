function newModel = changeBiomass(model,path,table,number,nonlipids,lipids,FAs)

% number:		Column number in textfile to get biomass coefficients from
% nonlipids:	Number of non-lipid components (typically 33), to be able to
%				separate what goes in biomass equation and what goes in lipid,
%				FA and AC equations.
% lipids:		Number of lipid species (typically 11), see above.
% FAs:			Number of fatty acids (typically >6), see above.
	if nargin<7
		FAs=10;
	end
	if nargin<6
		lipids=11;
	end
	if nargin<5
		nonlipids=33;
	end
    if nargin<4
        number=1;
    end
    if nargin<3
        table='no';
    end
    if nargin < 2
        path = 'D:\Users\eduardk.NET\Dropbox\Postdoc\Model\Yarrowia_lipolytica\Exp_biomassWTc.txt';
    end
    if nargin < 1
        model=model;
    end
    rxnstring = buildBiomass(path,table,number,nonlipids,lipids,FAs);
    rxnsToAdd.rxns={'xBIOMASS';'xLIPID';'xPOOL_FA_EM';'xPOOL_FA_EN';'xPOOL_FA_LP';'xPOOL_FA_MI';'xPOOL_FA_MM';'xPOOL_AC_EM';'xPOOL_AC_LP';'xPOOL_AC_MM'};
    ids    = getIndexes(model,rxnsToAdd.rxns,'rxns');
    rxnsToAdd.equations = rxnstring;
    rxnsToAdd.rxnNames  = model.rxnNames(ids);
    rxnsToAdd.lb        = [0;0;-1000;-1000;-1000;-1000;-1000;-1000;-1000;-1000];
    rxnsToAdd.ub        = [1000;1000;1000;1000;1000;1000;1000;1000;1000;1000];
    rxnsToAdd.c         = [0;0;0;0;0;0;0;0;0;0];
%      newModel=removeRxns(model,ids);
%      newModel=addRxns(newModel,rxnsToAdd,3);
     newModel=changeRxns(model,rxnsToAdd.rxns,rxnsToAdd.equations,3,{},true);
%      sol=solveLP(model,1)
%         newModel=changeRxns(model,rxnsToAdd.rxns(1),rxnsToAdd.equations(1),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(2),rxnsToAdd.equations(2),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(3),rxnsToAdd.equations(3),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(4),rxnsToAdd.equations(4),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(5),rxnsToAdd.equations(5),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(6),rxnsToAdd.equations(6),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(7),rxnsToAdd.equations(7),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(8),rxnsToAdd.equations(8),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(9),rxnsToAdd.equations(9),3);
%         sol_new=solveLP(newModel,1)
%         newModel=changeRxns(newModel,rxnsToAdd.rxns(10),rxnsToAdd.equations(10),3);
%         sol_new=solveLP(newModel,1)
        

end

function rxnstring = buildBiomass(path,table,number,nonlipids,lipids,FAs)
    %Reads a tab delimited file with met name - molecular weight - weight
    %fraction and stoichiometric coefficients and construct a BOF string
    if nargin<3
        number=1;
    end
    if nargin<2
        table='yes';
    end
    if nargin < 1
        path = 'D:\Users\eduardk.NET\Dropbox\Postdoc\Model\RandomSampling\Exp_biomass.txt';
    end
    fid         = fopen(path,'r');
    if strcmp(table,'yes');
        scan = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter','\t');
        BOFmets     = scan{1};
        BOFcoeff    = scan{number+1};
    elseif strcmp(table,'no');
        scan        = textscan(fid,'%s%f','Delimiter','\t');
        BOFmets     = scan{1};
        BOFcoeff    = scan{2};
    end
    rxnstring   = '';
    for i = 1:(nonlipids+1) %34
        rxnstring = horzcat(rxnstring,' + ',num2str(BOFcoeff(i)),' ',BOFmets{i});
    end
    ADPcoeff = BOFcoeff(strncmp('ATP[cy]',BOFmets,7));
    RHS = horzcat(' => biomass[cy] + ',num2str(ADPcoeff),' ADP[cy] + ',num2str(ADPcoeff),' phosphate[cy]');
    rxnstring = horzcat(rxnstring(4:end),RHS);
    
    tempstring='';
    for i = nonlipids+2:nonlipids+lipids+1%35:45
        tempstring = horzcat(tempstring,' + ',num2str(BOFcoeff(i)),' ',BOFmets{i});
    end
    rxnstring={rxnstring;horzcat(tempstring(4:end),' => 1000 lipids[cy]')};
    
    
    comps={'[em]','[en]','[lp]','[mi]','[mm]'};
    for j=1:length(comps)
        tempstring='';
        ratiocount=0;
    for i = nonlipids+lipids+FAs+2:nonlipids+lipids+FAs+FAs+1%58:69
        tempstring = horzcat(tempstring,' + ',num2str(BOFcoeff(i)),' ',BOFmets{i},comps{j});
    end
%    rxnstring{j+2}=horzcat(tempstring(4:end),' <=> fatty acid',comps{j});
    totalratio=sum(BOFcoeff(nonlipids+lipids+FAs+2:nonlipids+lipids+FAs+FAs+1));
    rxnstring{j+2}=horzcat(tempstring(4:end),' <=> ',num2str(totalratio),' fatty acid',comps{j});    
    end
    
    comps={'[em]','[lp]','[mm]'};
    for j=1:length(comps)
        tempstring='';
        ratiocount=0;
    for i = nonlipids+lipids+2:nonlipids+lipids+FAs+1%46:57
        tempstring = horzcat(tempstring,' + ',num2str(BOFcoeff(i)),' ',BOFmets{i},comps{j});
    end
%    rxnstring{j+7}=horzcat(tempstring(4:end),' <=> acyl-CoA',comps{j});
    totalratio=sum(BOFcoeff(nonlipids+lipids+2:nonlipids+lipids+FAs+1));
    rxnstring{j+7}=horzcat(tempstring(4:end),' <=> ',num2str(totalratio),' acyl-CoA',comps{j});
    end
end
