function out = discovery_symptom_associations(source)
% SOURCE is a private participant-level workbook; returns aggregate statistics only.
% Supplementary check of existing features; original inputs and manuscript unchanged.
% Primary: Spearman of OLS residuals on six manuscript covariates.
% P values use MATLAB corr on residuals, matching the prior validation procedure.

T=readtable(source,'VariableNamingRule','preserve');
roi=arrayfun(@(k)sprintf('MINDFIND_node_%03d',k),1:210,'UniformOutput',false);
X=T{:,roi};
D=table(string(T.sub),T.('组别（抑郁组=1，健康组=2）'),T.sex,T.age,...
 T.('独生子女（是=1，否=0）'),T.('单亲家庭（是=1.否=0）'),...
 T.mean_FD_Power_mm,T.TIV_mm3,T.SCL90,T.SDS,T.PHQ,...
 'VariableNames',{'ID','Group','sex','age','only','single','FD','TIV','SCL90','SDS','PHQ'});
covs={'sex','age','only','single','FD','TIV'};
keep=all(isfinite(X),2)&all(isfinite(D{:,['Group',covs]}),2);
D=D(keep,:);X=X(keep,:);
assert(height(D)==102 && sum(D.Group==1)==70 && sum(D.Group==2)==32);
assert(numel(unique(D.ID))==height(D));
assert(isequal(sort(unique(D.sex)),[1;2]));
assert(isequal(sort(unique(D.only)),[0;1]) && isequal(sort(unique(D.single)),[0;1]));
hc=find(D.Group==2);typ=nan(height(D),1);
for k=1:height(D)
 if D.Group(k)==2,reference=mean(X(setdiff(hc,k),:),1);else,reference=mean(X(hc,:),1);end
 typ(k)=corr(X(k,:)',reference');
end
D.typ=typ;
assert(all(isfinite(typ)) && all(abs(typ)<1));
out=struct('status','analysis of established cohort',...
 'n',height(D),'n_StD',70,'n_HC',32,...
 'template','32 discovery HC; leave-one-out for HC; fixed across symptom tests',...
 'covariates',{covs},'method','OLS residuals followed by two-sided Spearman; same procedure as validation',...
 'p_method','MATLAB corr(Type=Spearman) applied to residuals; not a rank-first partial Spearman test',...
 'FDR_family',{{'SCL90','SDS','PHQ'}},'inference','covariate-adjusted symptom associations in the established cohort');
scales={'SCL90','SDS','PHQ'};modes={'unadjusted','six_covariates'};
for mi=1:numel(modes)
 mode=modes{mi};ps=nan(3,1);rs=ps;results=struct;
 for si=1:numel(scales)
  scale=scales{si};A=D(isfinite(D.(scale)),:);xx=A.typ;yy=A.(scale);
  if mi>1
   cv=covs;
   C=A{:,cv};C=(C-mean(C))./std(C);design=[ones(height(A),1),C];
   assert(rank(design)==size(design,2));
   lm_x=fitlm(C,xx);lm_y=fitlm(C,yy);xx=lm_x.Residuals.Raw;yy=lm_y.Residuals.Raw;
   orthx=max(abs(design'*xx));orthy=max(abs(design'*yy));
   assert(orthx<1e-8 && orthy<1e-7);
   condition=cond(design);model_df=height(A)-size(design,2);
  else
   orthx=NaN;orthy=NaN;condition=NaN;model_df=NaN;
  end
  [rs(si),ps(si)]=corr(xx,yy,'Type','Spearman');
  results.(scale)=struct('n',height(A),'rho',rs(si),'p',ps(si),...
   'residual_model_df',model_df,'condition_number',condition,...
   'max_orthogonality_error_x',orthx,'max_orthogonality_error_y',orthy);
 end
 [sp,ix]=sort(ps);q=flipud(cummin(flipud(sp.*3./(1:3)')));adj=ps;adj(ix)=min(q,1);
 for si=1:numel(scales)
  results.(scales{si}).p_BH=adj(si);
  if si<=2,fprintf('%s %s n=%d rho=%+.9f p=%.9f BH=%.9f\n',mode,scales{si},results.(scales{si}).n,rs(si),ps(si),adj(si));end
 end
 out.(mode)=results;
end
assert(abs(out.unadjusted.SCL90.rho-(-0.24285916349020747))<1e-12);
assert(abs(out.unadjusted.SDS.rho-(-0.2328225440625537))<1e-12);

end
