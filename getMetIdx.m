function [ indexes ] = getMetIdx(model,met,comp,index)
%getMetIdx Obtains metabolite ID or index of given metabolite name and
%compartment
%   model       model structure (RAVEN)
%   met         string of metabolite name
%   comp        string of metabolite compartment, without square brackets
%   index       string indicating whether metabolite 'index' in
%               metabolite.mets should be given, or the actual metabolite
%               id. (optional, default: 'id')
%
%   2017-01-25      Eduard Kerkhoven (eduardk@chalmers.se)
if nargin<4
    index='id';
end
    metIdxs = find(strcmp(met,model.metNames));
    cmpidx = find(strcmp(comp,model.comps(model.metComps(metIdxs))));
    if strcmp(index,'index')
        indexes = metIdxs(cmpidx);
    else
        indexes = model.mets(metIdxs(cmpidx));
end
end

