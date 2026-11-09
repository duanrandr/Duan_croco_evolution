#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
# Step 1: Locate repository inputs and create a separate output directory.
code_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE))))
repo <- normalizePath(file.path(code_dir, '..', '..'))
out <- file.path(repo, 'results', 'reproduced', 'figure3')
dir.create(file.path(out, 'source_data'), recursive=TRUE, showWarnings=FALSE)
src <- file.path(repo, 'data', 'figure_source_data', 'figure3')
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
# Typeset superscripts/subscripts from Arial glyphs, avoiding Unicode fallback fonts.
ylabs <- expression(bold(paste('Specific growth rate (','d'^'-1',')')),
  bold(paste('N'['2'],' fixation / PON (','h'^'-1',')')),
  bold(paste('C fixation / POC (','h'^'-1',')')))
fes <- c('HFe','LFe')
fe_names <- c('Fe-replete  |  250 nM Fe','Fe-limited  |  5 nM Fe')
# Step 3: Define figure export functions.
export <- function(fig,stem,height) {
  w <- 180/25.4; h <- height/25.4
  ggsave(file.path(out,paste0(stem,'.pdf')),fig,device=vector_pdf,width=w,height=h,family='Arial',bg='white')
  ggsave(file.path(out,paste0(stem,'.svg')),fig,device=vector_svg,width=w,height=h,family='Arial',bg='white')
  if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,paste0(stem,'.png')),fig,device=ragg::agg_png,width=w,height=h,dpi=300,bg='white')
  if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,paste0(stem,'.tif')),fig,device=ragg::agg_tiff,width=w,height=h,dpi=600,compression='lzw',bg='white')
}
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
colour_guide <- guide_legend(order=1,nrow=1,override.aes=list(shape=15,alpha=1,size=2.4,stroke=.55,fill='white'))
base_prepare <- function(z) {
  z$Color <- as.character(z$History)
  z$center <- match(z$Temperature,c(23,28,32)) + c(-.26,0,.26)[match(z$History,c(23,28,32))]
  z
}
place_points <- function(z) {z$x <- z$center + c(-.055,0,.055)[z$Replicate];z}
# Main figure: selection history is encoded by color, measurement stage by shape.
d <- place_points(base_prepare(rd('Figure3_observations.csv')))
d$x <- d$center + c(-.04,0,.04)[d$Replicate]
s <- base_prepare(rd('Figure3_summary.csv'))
cld <- rd('letters.csv')
stopifnot(nrow(cld)==54, !anyDuplicated(paste(cld$Trait,cld$Treatment_Group)))
limits <- c(.58,.045,.092)
ticks <- list(seq(0,.5,.1),seq(0,.04,.01),seq(0,.08,.02))
# Step 4: Configure drawing styles, draw observations, and validate rendered layers.
plots <- list(); rendered_boxes <- list()
for(i in seq_along(traits)) for(j in seq_along(fes)) {
  dat <- d[d$Trait==traits[i]&d$Fe==fes[j],]; summ <- s[s$Trait==traits[i]&s$Fe==fes[j],]
  ann <- aggregate(Value~Treatment_Group+center,dat,max)
  ll <- cld[cld$Trait==traits[i],]
  ann$Letters <- ll$Letters[match(ann$Treatment_Group,ll$Treatment_Group)]
  ann$y <- ann$Value+.038*limits[i]
  stopifnot(nrow(ann)==9,all(!is.na(ann$Letters)),all(ann$y<limits[i]))
  q <- ggplot()+geom_vline(xintercept=c(1.5,2.5),colour='#808080',linewidth=.6,linetype='dashed')+
    geom_boxplot(data=dat,aes(x=center,y=Value,group=Treatment_Group,colour=Color,fill=Color),
      width=.20,linewidth=.34,alpha=.16,outlier.shape=NA,staplewidth=.45,show.legend=FALSE)+
    geom_point(data=dat,key_glyph=draw_key_temperature_or_role,aes(x=x,y=Value,colour=Color,shape=Role),size=1.9,stroke=.55,fill='white')+
    geom_text(data=ann,aes(x=center,y=y,label=Letters),family='Arial',fontface='bold',size=2.15,vjust=0,colour='#000000')+
    colour_scale+scale_fill_manual(values=temp_colors,guide='none')+
    scale_shape_manual(values=c('Short-term transfer'=16,'Long-term reference (no transfer)'=22),
      breaks=c('Short-term transfer','Long-term reference (no transfer)'),
      labels=c('Short-term transfer','Long-term selected (Tfinal)'),name=NULL)+
    scale_x_continuous(breaks=1:3,labels=c('23°C','28°C','32°C'),limits=c(.52,3.49),expand=expansion(0))+
    scale_y_continuous(limits=c(0,limits[i]),breaks=ticks[[i]],expand=expansion(mult=c(0,.003)),
      labels=if(i==1) function(v)sprintf('%.1f',v) else function(v)sprintf('%.2f',v))+
    guides(colour=colour_guide,shape=guide_legend(order=2,nrow=1,override.aes=list(colour='black',size=2.4,stroke=.55,fill='white')))+
    labs(y=if(j==1)ylabs[i] else NULL,x=if(i==3)'Measurement temperature' else NULL,
      title=if(i==1)fe_names[j] else NULL,tag=LETTERS[(i-1)*2+j])+common_theme
  built <- ggplot_build(q);bx <- built$data[[2]];pts <- built$data[[3]]
  stopifnot(nrow(pts)==nrow(dat), all(is.finite(pts$y)), all(abs(pts$y-dat$Value)<1e-12), nrow(built$data[[4]])==nrow(ann))
  for(g in unique(dat$Treatment_Group)) {
    gg <- dat[dat$Treatment_Group==g,]; row <- bx[abs(bx$x-gg$center[1])<1e-8,]
    stopifnot(nrow(row)==1, all(abs(gg$x-gg$center)<=.10))
    rendered_boxes[[length(rendered_boxes)+1]] <- data.frame(Treatment_Group=g,Trait=traits[i],Fe=fes[j],n=nrow(gg),
      Q1=row$lower,Median=row$middle,Q3=row$upper,Whisker_Lower=row$ymin,Whisker_Upper=row$ymax)
  }
  plots[[length(plots)+1]] <- q
}
write.csv(do.call(rbind,rendered_boxes),file.path(out,'source_data/Figure3_box_statistics.csv'),row.names=FALSE)
fig <- wrap_plots(G=guide_area(),A=plots[[1]],B=plots[[2]],C=plots[[3]],D=plots[[4]],E=plots[[5]],F=plots[[6]],
  design='GG\nAB\nCD\nEF',heights=c(.27,1.08,1,1.08),guides='collect')+
  plot_annotation(theme=theme(plot.margin=margin(4,4,3,3)))
export(fig,'Figure3_revised',162)
cat('Saved Figure 3 boxplots with original observations and verified significance letters.\n')
