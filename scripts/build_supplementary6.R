#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
# Step 1: Locate repository inputs and create a separate output directory.
code_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE))))
repo <- normalizePath(file.path(code_dir, '..', '..'))
out <- file.path(repo, 'results', 'reproduced', 'supplementary_6')
dir.create(file.path(out, 'source_data'), recursive=TRUE, showWarnings=FALSE)
src <- file.path(repo, 'data', 'figure_source_data', 'supplementary_6')
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
stars <- function(p) ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','ns')))
# Step 3: Define figure export functions.
export <- function(fig,stem,height) {
  w <- 180/25.4; h <- height/25.4
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
bracket_layers <- function(p,yspan) {
  br <- rbind(data.frame(x=p$x0,xend=p$x1,y=p$y,yend=p$y),
              data.frame(x=p$x0,xend=p$x0,y=p$y,yend=p$y-.015*yspan),
              data.frame(x=p$x1,xend=p$x1,y=p$y,yend=p$y-.015*yspan))
  list(geom_segment(data=br,aes(x=x,xend=xend,y=y,yend=yend),linewidth=.25,colour='#3c3c3c'),
       geom_text(data=p,aes(x=mid,y=label_y,label=label),family='Arial',fontface='bold',size=2.35,vjust=0,colour='#000000'))
}

# Validate the actual rendered box layer against type-7 quartiles and 1.5-IQR whiskers.
# All original observations are overlaid, including any beyond the box or whiskers.
audit_boxes <- function(plot,dat,box_layer,point_layer,width) {
  build<-ggplot_build(plot);bx<-build$data[[box_layer]];pt<-build$data[[point_layer]]
  stopifnot(nrow(pt)==nrow(dat),all(is.finite(pt$y)),all(abs(pt$y-dat$Value)<1e-12),all(abs(dat$x-dat$center)<=width/2))
  stopifnot(all(pt$size==1.9),all(pt$stroke==.55),all(pt$shape==16))
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

# Supplement: descriptive percentage distributions without significance annotations.

percent_ylabs <- expression(bold(atop('Growth rate','Change (%)')),
  bold(atop('N'['2']*' fixation / PON','Change (%)')),
  bold(atop('C fixation / POC','Change (%)')))
d <- rd('Source_temperature_change_observations.csv')
route_order <- c('23-28','23-32','28-32','32-28','32-23','28-23')
prep_pc <- function(z) {z$center<-match(paste(z$History,z$Temperature,sep='-'),route_order);z$Color<-as.character(z$History);z}
d<-place_points(prep_pc(d));plots<-list();pc_boxes<-list();pc_points<-0
d$x<-d$center+c(-.10,0,.10)[d$Replicate]
for(i in seq_along(traits)) for(j in seq_along(fes)) {
  dat<-d[d$Trait==traits[i]&d$Fe==fes[j],]
  all_d<-d[d$Trait==traits[i],]
  lo<-floor(min(all_d$Value,0)/25)*25-10
  hi<-ceiling(max(all_d$Value,0)/25)*25+25
  span<-hi-lo
  heading<-data.frame(x=c(2,5),y=rep(hi-.04*span,2),text=c('Warming','Cooling'))
  q<-ggplot()+geom_hline(yintercept=0,colour='#808080',linewidth=.6,linetype='dashed')+
    geom_vline(xintercept=3.5,colour='#808080',linewidth=.6,linetype='dashed')+
    geom_boxplot(data=dat,aes(x=center,y=Value,group=Treatment_Group,colour=Color,fill=Color),
      width=.52,linewidth=.34,alpha=.16,outlier.shape=NA,staplewidth=.45,show.legend=FALSE)+
    geom_point(data=dat,key_glyph=draw_key_temperature_or_role,aes(x=x,y=Value,colour=Color),shape=16,size=1.9,stroke=.55,fill='white')+
    geom_text(data=heading,aes(x=x,y=y,label=text),family='Arial',fontface='bold',size=2.3,vjust=1,colour='#000000')+
    colour_scale+scale_fill_manual(values=temp_colors,guide='none')+guides(colour=colour_guide)+
    scale_x_continuous(breaks=1:6,labels=c('23→28','23→32','28→32','32→28','32→23','28→23'),limits=c(.55,6.50),expand=expansion(0))+
    scale_y_continuous(limits=c(lo,hi),expand=expansion(0),breaks=function(x)pretty(x,n=4))+
    labs(y=if(j==1)percent_ylabs[i] else NULL,
         x=if(i==3)'Selection → measurement temperature (°C)' else NULL,
         title=if(i==1)fe_names[j] else NULL,tag=LETTERS[(i-1)*2+j])+common_theme+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=7.5,face='bold',colour='#000000'))
  # Only warming/cooling headings remain as text layers; no paired significance letters.
  stopifnot(length(q$layers)==5,nrow(ggplot_build(q)$data[[5]])==2)
  audit<-audit_boxes(q,dat,3,4,.52);pc_boxes[[length(pc_boxes)+1]]<-audit$boxes;pc_points<-pc_points+audit$points
  plots[[length(plots)+1]]<-q
}
fig<-wrap_plots(G=guide_area(),A=plots[[1]],B=plots[[2]],C=plots[[3]],D=plots[[4]],E=plots[[5]],F=plots[[6]],
  design='GG\nAB\nCD\nEF',heights=c(.16,1.08,1,1.14),guides='collect')+
  plot_annotation(theme=theme(plot.margin=margin(3,4,3,3)))
export(fig,'SuppFig6_source_temperature_change',152)
write.csv(do.call(rbind,pc_boxes),file.path(out,'source_data/Source_temperature_change_box_statistics.csv'),row.names=FALSE)
stopifnot(pc_points==107)

cat('Saved unchanged descriptive source-temperature percentage figure.\n')
