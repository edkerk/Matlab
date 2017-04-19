function prntflx(modelName,solName,outputFile)
% prntflx
%   Print flux distribution from RAVEN simulation to text file
%
%   modelName     name of simulated model (default=model)
%   solName       name of simulation solution vector (default=sol.x)
%   outputFile    name of output file (default=tmp.tab)



if nargin<3
    outputFile='tmp.tab';
end

if nargin<2
    solName=sol.x;
end

if nargin<1
    modelname=model;
end

printFluxes(modelName,solName,false,[],outputFile,'%rxnID\t%rxnName\t%eqn\t%flux\n');