#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
# Step 1: Locate repository inputs and create a separate output directory.
code_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE))))
repo <- normalizePath(file.path(code_dir, '..', '..'))
out <- file.path(repo, 'results', 'reproduced', 'supplementary_7')
dir.create(file.path(out, 'source_data'), recursive=TRUE, showWarnings=FALSE)
src <- file.path(repo, 'data', 'figure_source_data', 'supplementary_7')
# Step 2: Load export devices and read the supplied figure inputs.
source(file.path(code_dir, 'vector_export.R'))
# Temperature colors use rectangular swatches; role keys retain their point symbols.
draw_key_temperature_or_role <- function(data, params, size) {
  if (!is.null(data$shape) && !is.na(data$shape) && data$shape == 15) {
    grid::rectGrob(width=grid::unit(3.2,'mm'),height=grid::unit(1.5,'mm'),
      gp=grid::gpar(fill=data$colour,col=NA))
  } else {
    ggplot2::draw_key_point(data, params, size)
  }
}
rd <- function(x) read.csv(file.path(src,x),stringsAsFactors=FALSE)
temp_colors <- c('23'='#4a1f1f','28'='#115cc8','32'='#ff1a1a','T0'='#676767')
traits <- c('Growth_rates','Nfix_PON','Cfix_POC')
ylabs <- c('Specific growth rate (d⁻¹)','N₂ fixation / PON (h⁻¹)','C fixation / POC (h⁻¹)')
fes <- c('HFe','LFe')
fe_names <- c('Fe-replete  |  250 nM Fe','Fe-limited  |  5 nM Fe')
# Step 3: Define figure export functions.
export <- function(fig,stem,height) {
  w <- 120/25.4; h <- height/25.4
  ggsave(file.path(out,paste0(stem,'.pdf')),fig,device=vector_pdf,width=w,height=h,family='Arial',bg='white')
  ggsave(file.path(out,paste0(stem,'.svg')),fig,device=vector_svg,width=w,height=h,family='Arial',bg='white')
  if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,paste0(stem,'.png')),fig,device=ragg::agg_png,width=w,height=h,dpi=300,bg='white')
  if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,paste0(stem,'.tif')),fig,device=ragg::agg_tiff,width=w,height=h,dpi=600,compression='lzw',bg='white')
}
# Step 4: Configure drawing styles, draw observations, and validate rendered layers.
common_theme <- theme_classic(base_family='Arial',base_size=7) +
  theme(text=element_text(family='Arial',face='bold',colour='#000000'),
        axis.line=element_line(colour='#000000',linewidth=.60),axis.ticks=element_line(colour='#000000',linewidth=.50),
        axis.ticks.length=grid::unit(1,'mm'),axis.text=element_text(size=7.5,colour='#000000',face='bold'),
        axis.title.y=element_text(size=8.5,face='bold',colour='#000000',margin=margin(r=5)),axis.title.x=element_text(size=8.5,face='bold',colour='#000000',margin=margin(t=4)),
        plot.title=element_text(size=8.5,face='bold',colour='#000000',hjust=.5,margin=margin(b=4)),
        plot.tag=element_text(size=9,face='bold',colour='#000000'),plot.tag.position='topleft',
        plot.margin=margin(3,5,2,3),legend.position='top',legend.box='vertical',
        legend.title=element_text(size=7.5,face='bold',colour='#000000'),legend.text=element_text(size=7.5,face='bold',colour='#000000'),
        legend.key.width=grid::unit(4,'mm'),legend.key.height=grid::unit(3.2,'mm'),
        legend.margin=margin(0,0,0,0),legend.spacing.y=grid::unit(0,'mm'))
colour_scale <- scale_color_manual(values=temp_colors,breaks=c('23','28','32'),labels=c('23°C','28°C','32°C'),name='Long-term selection temperature')
colour_guide <- guide_legend(order=1,nrow=1,override.aes=list(shape=15,alpha=1,size=2.4,stroke=.55))
base_prepare <- function(z) {
  z$Color <- as.character(z$History)
  z$center <- match(z$Temperature,c(23,28,32)) + c(-.26,0,.26)[match(z$History,c(23,28,32))]
  z
}
place_points <- function(z) {z$x <- z$center + c(-.055,0,.055)[z$Replicate];z}
# Validate the actual rendered box layer against type-7 quartiles and 1.5-IQR whiskers.
# All original observations are overlaid, including any beyond the box or whiskers.
audit_boxes <- function(plot,dat,box_layer,point_layer,width) {
  build<-ggplot_build(plot);bx<-build$data[[box_layer]];pt<-build$data[[point_layer]]
  stopifnot(nrow(pt)==nrow(dat),all(is.finite(pt$y)),all(abs(pt$y-dat$Value)<1e-12),all(abs(dat$x-dat$center)<=width/2))
  stopifnot(all(pt$size==1.9),all(pt$stroke==.55))
  for(k in seq_along(plot$layers)) {
    if(inherits(plot$layers[[k]]$geom,'GeomHline') || inherits(plot$layers[[k]]$geom,'GeomVline')) {
      ln<-build$data[[k]];stopifnot(all(ln$colour=='#808080'),all(ln$linewidth==.6),all(ln$linetype=='dashed'))
    }
  }
  rows<-lapply(unique(dat$Treatment_Group),function(g) {
    z<-dat[dat$Treatment_Group==g,];b<-bx[abs(bx$x-z$center[1])<1e-10,]
    qu<-quantile(z$Value,c(.25,.5,.75),type=7);iq<-qu[3]-qu[1]
    wh<-range(z$Value[z$Value>=qu[1]-1.5*iq & z$Value<=qu[3]+1.5*iq])
    stopifnot(nrow(b)==1,max(abs(c(b$lower,b$middle,b$upper)-qu))<1e-10,
      max(abs(c(b$ymin,b$ymax)-wh))<1e-10,all(z$x>=b$xmin),all(z$x<=b$xmax))
    data.frame(Treatment_Group=g,Trait=z$Trait[1],Fe=z$Fe[1],n=nrow(z),
      Q1=b$lower,Median=b$middle,Q3=b$upper,Whisker_Lower=b$ymin,Whisker_Upper=b$ymax)
  })
  list(boxes=do.call(rbind,rows),points=nrow(pt))
}


# Matched-condition comparison: true ancestral T0 and transfers measured at 28 C HFe.
# The original box/point/color/font style is retained. Only the comparison scope
# and its statistical family change from the original all-condition figure.
d<-rd('T0_growth_observations.csv');cld<-rd('T0_growth_letters.csv')
cld<-cld[cld$Trait=='Growth_rates',]
groups<-c('T0_28HFe','23-28HFe','32-28HFe')
stopifnot(nrow(cld)==3, setequal(cld$Treatment_Group,groups),
  !anyDuplicated(cld$Treatment_Group),all(nchar(cld$Letters)>0),
  nrow(d)==9,setequal(d$Treatment_Group,groups),all(d$Fe=='HFe'),all(d$Temperature==28))
d$Color<-ifelse(d$Treatment_Group=='T0_28HFe','T0',as.character(d$History))
d$center<-match(d$Treatment_Group,groups)
d$x<-d$center+c(-.055,0,.055)[d$Replicate]
lim<-.58
lab<-aggregate(Value~Treatment_Group+Fe+center,d,max)
lab$Letters<-cld$Letters[match(lab$Treatment_Group,cld$Treatment_Group)]
lab$label_y<-lab$Value+.035*lim
q<-ggplot()+
  geom_boxplot(data=d,aes(x=center,y=Value,group=Treatment_Group,colour=Color,fill=Color),
    width=.28,linewidth=.34,alpha=.16,outlier.shape=NA,staplewidth=.45,show.legend=FALSE)+
  geom_point(data=d,aes(x=x,y=Value,colour=Color,shape=Role),size=1.9,stroke=.55,fill='white',show.legend=FALSE)+
  geom_text(data=lab,aes(x=center,y=label_y,label=Letters),family='Arial',fontface='bold',size=2.15,vjust=0,colour='#000000')+
  scale_color_manual(values=temp_colors,guide='none')+scale_fill_manual(values=temp_colors,guide='none')+
  scale_shape_manual(values=c('Ancestral (T0)'=23,'Short-term transfer'=16),guide='none')+
  scale_x_continuous(breaks=1:3,labels=c('Ancestral (T0)','23→28°C','32→28°C'),limits=c(.55,3.45),expand=expansion(0))+
  scale_y_continuous(limits=c(0,lim),breaks=seq(0,.5,.1),expand=expansion(0),labels=function(v)sprintf('%.1f',v))+
  labs(y=expression(bold(paste('Specific growth rate (','d'^'-1',')'))),
       x='All cultures measured at 28°C',title='Fe-replete  |  250 nM Fe')+common_theme+
  theme(plot.margin=margin(7,10,7,7))
audit<-audit_boxes(q,d,1,2,.28)
stopifnot(audit$points==9,nrow(ggplot_build(q)$data[[3]])==3)
export(q,'SuppFig7_T0_growth_comparison',80)
write.csv(audit$boxes,file.path(out,'source_data/T0_growth_box_statistics.csv'),row.names=FALSE)
write.csv(lab,file.path(out,'source_data/T0_growth_letter_positions.csv'),row.names=FALSE)
cat('Verified and saved three matched-condition growth groups, nine observations, and new three-group letters.\n')
