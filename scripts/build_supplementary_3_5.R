#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
# Step 1: Locate repository inputs and create a separate output directory.
code_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE))))
repo <- normalizePath(file.path(code_dir, '..', '..'))
out <- file.path(repo, 'results', 'reproduced', 'supplementary_3_5')
dir.create(file.path(out, 'source_data'), recursive=TRUE, showWarnings=FALSE)
src <- file.path(repo, 'data', 'figure_source_data', 'supplementary_3_5')
# Step 2: Load export devices and read the supplied figure inputs.
source(file.path(code_dir, 'vector_export.R'))
# Export each supplementary figure in its own graphics session.
trait_choice<-Sys.getenv('CROCO_SUPP_TRAIT')
if(!nzchar(trait_choice)) {
  this_script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)))
  for(k in 1:3) {
    status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(this_script),env=paste0('CROCO_SUPP_TRAIT=',k))
    stopifnot(status==0)
  }
  boxes<-lapply(1:3,function(k) read.csv(file.path(out,paste0('source_data/Supplementary_box_statistics_',k,'.csv'))))
  write.csv(do.call(rbind,boxes),file.path(out,'source_data/Supplementary_box_statistics.csv'),row.names=FALSE)
  cat('Saved all three descriptive supplementary figures in the requested formats.\n')
  quit(save='no')
}
stopifnot(trait_choice %in% c('1','2','3'))
# Temperature colors use rectangular swatches; role keys retain their point symbols.
draw_key_temperature_or_role <- function(data, params, size) {
  if (!is.null(data$shape) && !is.na(data$shape) && data$shape == 15) {
    grid::rectGrob(width=grid::unit(3.2,'mm'),height=grid::unit(1.5,'mm'),
      gp=grid::gpar(fill=data$colour,col=NA))
  } else {
    ggplot2::draw_key_point(data, params, size)
  }
}
d<-read.csv(file.path(src,'All_treatment_observations.csv'),stringsAsFactors=FALSE)
stopifnot(length(unique(d$Treatment_Group))==25, all(d$Treatment_Group[grepl('^T0_',d$Treatment_Group)]=='T0_28HFe'))
colors<-c('T0'='#676767','23'='#4a1f1f','28'='#115cc8','32'='#ff1a1a')
roles<-c('Ancestral (T0)','Initial acclimation','Long-term selected','Short-term transfer')
d$Role<-factor(d$Role,levels=roles)
d$Color<-ifelse(d$Role==roles[1],'T0',as.character(d$History))
# Keep the measurement-temperature layout, adding one separate IA position.
# Five positions are needed only at 28 HFe (ancestral + IA + three selected histories).
position_key<-ifelse(d$Role==roles[1],'T0',ifelse(d$Role==roles[2],'IA',as.character(d$History)))
offset<-c(-.27,-.09,.09,.27)[match(position_key,c('IA','23','28','32'))]
at_ancestor<-d$Temperature==28 & d$Fe=='HFe'
offset[at_ancestor]<-c(-.36,-.18,0,.18,.36)[match(position_key[at_ancestor],c('T0','IA','23','28','32'))]
d$center<-match(d$Temperature,c(23,28,32))+offset
d$x<-d$center+c(-.04,0,.04)[d$Replicate]
traits<-c('Growth_rates','Nfix_PON','Cfix_POC')
stems<-c('SuppFig3_growth_responses','SuppFig4_nfix_responses','SuppFig5_cfix_responses')
ylabs<-expression(bold(paste('Specific growth rate (','d'^'-1',')')),
  bold(paste('N'['2'],' fixation / PON (','h'^'-1',')')),
  bold(paste('C fixation / POC (','h'^'-1',')')))
limits<-c(.52,.042,.083)
ticks<-list(seq(0,.5,.1),seq(0,.04,.01),seq(0,.08,.02))
fes<-c('HFe','LFe');fe_names<-c('Fe-replete  |  250 nM Fe','Fe-limited  |  5 nM Fe')
# Step 3: Draw observations and validate the rendered layers.
rendered_boxes<-list()
for(i in as.integer(trait_choice)) {
  panels<-list()
  for(j in seq_along(fes)) {
    dat<-d[d$Trait==traits[i]&d$Fe==fes[j],]
    q<-ggplot()+geom_vline(xintercept=c(1.5,2.5),linewidth=.6,colour='#808080',linetype='dashed')+
      geom_boxplot(data=dat,aes(x=center,y=Value,group=Treatment_Group,colour=Color,fill=Color),
        width=.15,linewidth=.33,alpha=.16,staplewidth=.45,outlier.shape=NA,show.legend=FALSE)+
      geom_point(data=dat,key_glyph=draw_key_temperature_or_role,aes(x=x,y=Value,colour=Color,shape=Role),size=1.9,stroke=.55,fill='white',show.legend=TRUE)+
      scale_color_manual(values=colors,breaks=c('23','28','32'),labels=c('23°C','28°C','32°C'),name='Initial / selection temperature')+
      scale_fill_manual(values=colors,guide='none')+
      scale_shape_manual(values=c('Ancestral (T0)'=23,'Initial acclimation'=24,'Long-term selected'=22,'Short-term transfer'=16),breaks=roles,labels=c('Ancestral (T0)','Initial acclimation','Long-term selected (Tfinal)','Short-term transfer'),name=NULL,drop=FALSE)+
      scale_x_continuous(breaks=1:3,labels=c('23°C','28°C','32°C'),limits=c(.48,3.54),expand=expansion(0))+
      scale_y_continuous(limits=c(0,limits[i]),breaks=ticks[[i]],expand=expansion(mult=c(0,.003)),
        labels=if(i==1) function(v)sprintf('%.1f',v) else function(v)sprintf('%.2f',v))+
      guides(colour=guide_legend(order=1,nrow=1,override.aes=list(shape=15,alpha=1,size=2.4,stroke=.55)),
        shape=guide_legend(order=2,nrow=1,override.aes=list(colour=c('#676767','black','black','black'),size=2.4,stroke=.55,fill='white')))+
      labs(x='Measurement temperature',y=if(j==1)ylabs[i] else NULL,title=fe_names[j],tag=LETTERS[j])+
      theme_classic(base_size=7,base_family='Arial')+
      theme(axis.line=element_line(linewidth=.60,colour='#000000'),axis.ticks=element_line(linewidth=.50,colour='#000000'),
        axis.ticks.length=grid::unit(1,'mm'),axis.text=element_text(size=7.5,face='bold',colour='#000000'),
        axis.title.y=element_text(size=8.5,face='bold',colour='#000000',margin=margin(r=5)),axis.title.x=element_text(size=8.5,face='bold',colour='#000000',margin=margin(t=4)),
        plot.title=element_text(size=8.5,face='bold',colour='#000000',hjust=.5,margin=margin(b=4)),
        plot.tag=element_text(size=9,face='bold',colour='#000000'),plot.tag.position='topleft',
        plot.margin=margin(3,5,2,3),legend.position='top',legend.box='vertical',
        legend.title=element_text(size=7.5,face='bold',colour='#000000'),legend.text=element_text(size=7.5,face='bold',colour='#000000'),
        legend.key.width=grid::unit(3.6,'mm'),legend.key.height=grid::unit(3.2,'mm'),
        legend.margin=margin(0,0,0,0),legend.spacing.y=grid::unit(0,'mm'))
    built<-ggplot_build(q);bx<-built$data[[2]];pts<-built$data[[3]]
    stopifnot(nrow(pts)==nrow(dat),all(is.finite(pts$y)),
      all(abs(pts$x-dat$x)<1e-12),all(abs(pts$y-dat$Value)<1e-12),
      length(q$layers)==3,!any(vapply(q$layers,function(z) inherits(z$geom,'GeomText'),logical(1))))
    for(g in unique(dat$Treatment_Group)) {
      gg<-dat[dat$Treatment_Group==g,];row<-bx[abs(bx$x-gg$center[1])<1e-8,]
      stopifnot(nrow(row)==1,all(abs(gg$x-gg$center)<=.15/2))
      rendered_boxes[[length(rendered_boxes)+1]]<-data.frame(Figure=stems[i],Trait=traits[i],Treatment_Group=g,Fe=fes[j],n=nrow(gg),
        Q1=row$lower,Median=row$middle,Q3=row$upper,Whisker_Lower=row$ymin,Whisker_Upper=row$ymax)
    }
    panels[[j]]<-q
  }
  fig<-wrap_plots(G=guide_area(),A=panels[[1]],B=panels[[2]],design='GG\nAB',heights=c(.27,1),guides='collect')+
    plot_annotation(theme=theme(plot.margin=margin(3,4,3,3)))
  w<-180/25.4;h<-79/25.4
# Step 4: Save recreated figures and numerical drawing checks.
  ggsave(file.path(out,paste0(stems[i],'.pdf')),fig,device=vector_pdf,width=w,height=h,family='Arial',bg='white')
  ggsave(file.path(out,paste0(stems[i],'.svg')),fig,device=vector_svg,width=w,height=h,family='Arial',bg='white')
  if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,paste0(stems[i],'.png')),fig,device=ragg::agg_png,width=w,height=h,dpi=300,bg='white')
  if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,paste0(stems[i],'.tif')),fig,device=ragg::agg_tiff,width=w,height=h,dpi=600,compression='lzw',bg='white')
}
write.csv(do.call(rbind,rendered_boxes),file.path(out,paste0('source_data/Supplementary_box_statistics_',trait_choice,'.csv')),row.names=FALSE)
cat('Saved descriptive Supplementary Figures 3, 4 and 5 with boxplots and original observations, without significance letters.\n')
