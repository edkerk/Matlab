function [ output_args ] = randSample_and_corr( model_ref,model_sample,RNA,protein )
%randSample_and_corr Perform random sampling and do correlation analysis
%with RNA and/or protein
%   Input:
%       model_ref           model structure, with contraints fixed for
%                           reference conditions
%       model_control       model structure, with contraints fixed for
%                           control conditions

% 
if ~exist('protein', 'var') ||  isempty(protein)
    protein = false;
end
if ~exist('RNA', 'var') ||  isempty(RNA)
    RNA = false;
end
    
% Perform random sampling on both models and set very low flux to zero and
% calculate z-scores
rs.ref=full(transpose(randomSampling(model_ref)));
rs.sample=full(transpose(randomSampling(model_sample)));
rs.ref(abs(rs.ref)<1e-10)=0;
rs.sample(abs(rs.sample)<1e-10)=0;
flux.mean.ref=mean(rs.ref);
flux.mean.sample=mean(rs.sample);
flux.sd.ref=std(rs.ref);
flux.sd.sample=std(rs.sample);
z=transpose((flux.mean.sample-flux.mean.ref)./(sqrt(var(rs.sample)+var(rs.ref))));


out=
% Prepare output structure


end

