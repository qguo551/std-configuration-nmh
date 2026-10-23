function out = validation_symptom_associations(source)
% SOURCE is a private participant-level workbook; returns aggregate statistics only.
% Diagnostic comparison on the established validation cohort and feature.
T=readtable(source,'VariableNamingRule','preserve');
roi=arrayfun(@(k)sprintf('MINDFIND_node_%03d',k),1:210,'UniformOutput',false);
X=T{:,roi};g=T.('组别（抑郁组=1，健康组=2）');
keep=all(isfinite(X),2)&isfinite(T.mean_FD_Power_mm)&isfinite(T.TIV_mm3)&isfinite(g);
T=T(keep,:);X=X(keep,:);g=g(keep);
assert(sum(g==1)==29 && sum(g==2)==16);
reference=mean(X(g==2,:),1);xd=X(g==1,:);typ=nan(29,1);
for k=1:29,typ(k)=corr(xd(k,:)',reference');end
T=T(g==1,:);
D=table(string(T.sub),T.sex,T.age,T.('独生子女（是=1，否=0）'),...
 T.('单亲家庭（是=1.否=0）'),T.mean_FD_Power_mm,T.TIV_mm3,T.SCL90,T.SDS,typ,...
 'VariableNames',{'ID','sex','age','only','single','FD','TIV','SCL90','SDS','typ'});
assert(numel(unique(D.ID))==29);
covs={'sex','age','only','single','FD','TIV'};
assert(all(isfinite(D{:,covs}),'all') && all(isfinite(D{:,{'typ','SCL90','SDS'}}),'all'));
assert(isequal(sort(unique(D.sex)),[1;2]));
assert(isequal(sort(unique(D.only)),[0;1]) && isequal(sort(unique(D.single)),[0;1]));
out=struct('n',29,'n_reference_HC',16,...
 'covariates',{covs},'method','Separate OLS residuals followed by two-sided Spearman, matching existing validation',...
 'p_method','MATLAB corr(Type=Spearman) on residuals; not rank-first partialcorr',...
 'multiplicity','Unadjusted p values, matching existing validation reporting',...
 'role','Post hoc diagnostic comparison; no manuscript update; no new exclusions');
for mode={'FD_TIV'}
 name=mode{1};cv=covs;if strcmp(name,'FD_TIV'),cv={'FD','TIV'};end
 if strcmp(name,'four_covariates'),cv={'sex','age','FD','TIV'};end
 if strcmp(name,'five_without_only'),cv={'sex','age','single','FD','TIV'};end
 if strcmp(name,'five_without_single'),cv={'sex','age','only','FD','TIV'};end
 out.model_covariates.(name)=cv;
 C=D{:,cv};C=(C-mean(C))./std(C);A=[ones(29,1),C];assert(rank(A)==size(A,2));
 for scale={'SCL90','SDS'}
  s=scale{1};lx=fitlm(C,D.typ);ly=fitlm(C,D.(s));rx=lx.Residuals.Raw;ry=ly.Residuals.Raw;
  assert(max(abs(A'*rx))<1e-8 && max(abs(A'*ry))<1e-7);
  [rho,p]=corr(rx,ry,'Type','Spearman');
  out.(name).(s)=struct('n',29,'rho',rho,'p',p,'residual_model_df',29-size(A,2),'design_rank',rank(A),'condition_number',cond(A));
  fprintf('%s %s n=29 rho=%+.12f p=%.12f\n',name,s,rho,p);
 end
end
assert(abs(out.FD_TIV.SCL90.rho-(-0.42413793103448283))<1e-12);
assert(abs(out.FD_TIV.SCL90.p-0.022711531070421634)<1e-12);
assert(abs(out.FD_TIV.SDS.rho-(-0.096059113300492618))<1e-12);

end
