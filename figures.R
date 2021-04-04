filenames=dir(pattern = '^simu_')


datalist_sa0.1=list('nmod'=NULL,'amod'=NULL)

method=c('MR-Egger','LDA MR-Egger','TWMR','RAPS','MR-LDP','LDMR','PLDMRa','PLDMR','PLDMRt')

for (i in 1:length(filenames))
{
  para=as.numeric(unlist(regmatches(filenames[i],regexec("simu_beta([0-9]+\\.+[0-9]+|[0-9]+)a([-]?[0-9]+\\.+[0-9]+|[0-9]+)sa([0-9]+\\.+[0-9]+|[0-9]+)LD([0-9]+\\.+[0-9]+|[0-9]+)n([0-9]+)p([0-9]+)\\.txt", filenames[i]))))
  beta=para[2]
  type=para[3]
  sa=para[4]
  LD=para[5]
  n=para[6]
  m=para[7]
    data=read.table(filenames[i],sep='\t',row.names=NULL,header=T)
    res_beta=unname(colMeans(data[,1:9]))
    res_se=unname(sqrt(colMeans((data[,1:9]-res_beta)^2)))
    res_lower=res_beta-res_se
    res_upper=res_beta+res_se
    res_min=unname(apply(data[1:9], 2, min))
    res_max=unname(apply(data[1:9], 2, max))
    res_pval=unname(colMeans(data[,10:18]<0.05))
    res=data.frame(Method=method,N=n,M=m,LD=paste('rho[g]==',LD,sep=''),Type=paste('mu[alpha]==',type,sep=''),beta_mean=res_beta,beta_se=res_se,beta_lower=res_lower,beta_upper=res_upper,
                   beta_min=res_min,beta_max=res_max,Power=res_pval)
    if (m>25) {
      if (beta==0) {
        datalist_sa0.1$nmod=rbind(res,datalist_sa0.1$nmod)
      } else {
        datalist_sa0.1$amod=rbind(res,datalist_sa0.1$amod)
      }
    }
}

library(ggplot2)
res=datalist_sa0.1$nmod
res=res[res$N>=1000&res$Method!='TWMR',]
res$M=factor(res$M,levels=c('50','75','100','200'))
res$LD=factor(res$LD,levels=c('rho[g]==0','rho[g]==0.3','rho[g]==0.6'))
res$Type=factor(res$Type,levels=c('mu[alpha]==0','mu[alpha]==-0.1','mu[alpha]==0.1'))
res$Method=factor(res$Method,levels=c('MR-LDP','RAPS','MR-Egger','LDA MR-Egger','LDMR','PLDMRa','PLDMR','PLDMRt'))
pd <- position_dodge(width=0.7)
ggplot(data=res,aes(x=M,y=Power,fill=Method))+geom_bar(stat='identity',position=position_dodge())+
  facet_grid(LD~Type,labeller=label_parsed)+scale_x_discrete(breaks=c('50','75','100','200'),labels=c('50','75','100','200'))+
  geom_hline(yintercept=0.05,col='red',lty=2)+geom_hline(yintercept=0.1,col='blue',lty=3)+xlab('Numbers of IVs')+ylab(label='Type I Error Rates')+theme(text = element_text(size=10))
ggsave('T1em.png',device='png',dpi='print',width=18,height=12,units='cm')

ggplot(data=res, aes(x=M, y=beta_mean,color=Method))+geom_errorbar(aes(ymin=beta_lower,ymax=beta_upper),position=pd)+geom_point(position=pd)+
  facet_grid(LD~Type,labeller=label_parsed)+scale_x_discrete(breaks=c('50','75','100','200'),labels=c('50','75','100','200'))+
  ylab(label=expression(hat(beta)%+-%~SE))+xlab('Numbers of IVs')+theme(text = element_text(size=10))
ggsave('beta0m.eps',device='eps',dpi='print',width=18,height=12,units='cm')

res=datalist_sa0.1$amod
res=res[res$N>=1000,]
res=res[res$N>=1000&res$Method!='TWMR',]
res$M=factor(res$M,levels=c('50','75','100','200'))
res$LD=factor(res$LD,levels=c('rho[g]==0','rho[g]==0.3','rho[g]==0.6'))
res$Type=factor(res$Type,levels=c('mu[alpha]==0','mu[alpha]==-0.1','mu[alpha]==0.1'))
res$Method=factor(res$Method,levels=c('MR-LDP','RAPS','MR-Egger','LDA MR-Egger','LDMR','PLDMRa','PLDMR','PLDMRt'))
pd <- position_dodge(width=0.7)
ggplot(data=res,aes(x=M,y=Power,fill=Method))+geom_bar(stat='identity',position=position_dodge())+
  facet_grid(LD~Type,labeller=label_parsed)+scale_x_discrete(breaks=c('50','75','100','200'),labels=c('50','75','100','200'))+
  ylab(label='Powers')+xlab('Numbers of IVs')+theme(text = element_text(size=10))
ggsave('Powerm.png',device='png',dpi='print',width=18,height=12,units='cm')

ggplot(data=res, aes(x=M, y=beta_mean,color=Method))+geom_errorbar(aes(ymin=beta_lower,ymax=beta_upper),position=pd)+geom_point(position=pd)+
  facet_grid(LD~Type,labeller=label_parsed)+scale_x_discrete(breaks=c('50','75','100','200'),labels=c('50','75','100','200'))+
  ylab(label=expression(hat(beta)%+-%~SE))+xlab('Numbers of IVs')+theme(text = element_text(size=10))
ggsave('beta0.05m.eps',device='eps',dpi='print',width=18,height=12,units='cm')

