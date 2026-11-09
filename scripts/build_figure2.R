#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
# Step 1: Locate repository inputs and create a separate output directory.
code_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE))))
repo <- normalizePath(file.path(code_dir, '..', '..'))
out <- file.path(repo, 'results', 'reproduced', 'figure2')
dir.create(file.path(out, 'source_data'), recursive=TRUE, showWarnings=FALSE)
src <- file.path(repo, 'data', 'figure_source_data', 'figure2')
# Step 2: Load export devices and read the supplied figure inputs.
source(file.path(code_dir, 'vector_export.R'))
d <- read.csv(file.path(src, 'Figure2_observations.csv'))
s <- read.csv(file.path(src, 'Figure2_summary.csv'))
cld <- read.csv(file.path(src, 'letters.csv'), stringsAsFactors=FALSE)
cld <- cld[c('Trait','Treatment_Group','Letters')]
stopifnot(nrow(cld)==39, all(c('Trait','Treatment_Group','Letters') %in% names(cld)),
          !anyDuplicated(cld[c('Trait','Treatment_Group')]),
          all(!is.na(cld$Letters) & nzchar(cld$Letters)))
stopifnot(length(unique(d$Treatment_Group))==13, all(d$Treatment_Group[grepl('^T0_',d$Treatment_Group)]=='T0_28HFe'))
stages <- c('Ancestral (T0)', 'Initial acclimation', 'Long-term selected')
temp_colors <- c('23'='#4a1f1f', '28'='#115cc8', '32'='#ff1a1a', 'T0'='#676767')
prepare <- function(z) {
  z$Stage <- factor(z$Stage, levels=stages)
  z$Fe <- factor(z$Fe, levels=c('HFe','LFe'))
  z$Temp <- factor(z$Temperature, levels=c(23,28,32))
  z$Color <- ifelse(z$Stage==stages[1],'T0',as.character(z$Temperature))
  # Retain the established temperature locations; separate all three stages at 28 HFe.
  offset <- ifelse(z$Stage==stages[2],-.19,.19)
  at_ancestor <- z$Temperature==28 & z$Fe=='HFe'
  offset[at_ancestor] <- c(-.29,0,.29)[match(z$Stage[at_ancestor],stages)]
  z$center <- match(z$Temperature,c(23,28,32)) + offset
  z
}
d <- prepare(d); s <- prepare(s)
# Points share their own box center; only horizontal displacement is applied.
d$x <- d$center + c(-.05,0,.05)[d$Replicate]
traits <- c('Growth_rates','Nfix_PON','Cfix_POC')
fes <- c('HFe','LFe')
fe_names <- c('Fe-replete  |  250 nM Fe','Fe-limited  |  5 nM Fe')
# Typeset superscripts/subscripts from Arial glyphs, avoiding Unicode fallback fonts.
ylabs <- expression(bold(paste('Specific growth rate (','d'^'-1',')')),
  bold(paste('N'['2'],' fixation / PON (','h'^'-1',')')),
  bold(paste('C fixation / POC (','h'^'-1',')')))
limits <- c(.55,.040,.080)
breaks <- list(seq(0,.5,.1), seq(0,.04,.01), seq(0,.08,.02))
# Step 3: Draw observations and validate the rendered layers.
plots <- list(); rendered_boxes <- list(); rendered_letters <- list()
for(i in seq_along(traits)) for(j in seq_along(fes)) {
  trait <- traits[i]
  dat <- d[d$Trait==trait & d$Fe==fes[j],]
  summ <- s[s$Trait==trait & s$Fe==fes[j],]
  annotation <- merge(summ,cld[cld$Trait==trait,],by=c('Trait','Treatment_Group'),sort=FALSE)
  stopifnot(nrow(annotation)==if(j==1)7 else 6,
            setequal(annotation$Treatment_Group,dat$Treatment_Group))
  annotation$y <- vapply(annotation$Treatment_Group,function(g) max(dat$Value[dat$Treatment_Group==g]),numeric(1)) + .035*limits[i]
  stopifnot(all(annotation$y < limits[i]))
  q <- ggplot() +
    geom_vline(xintercept=c(1.5,2.5),colour='#808080',linewidth=.6,linetype='dashed') +
    geom_boxplot(data=dat,aes(x=center,y=Value,group=Treatment_Group,colour=Color,fill=Color),
      width=.23,linewidth=.36,alpha=.16,outlier.shape=NA,staplewidth=.45,show.legend=FALSE) +
    geom_point(data=dat,aes(x=x,y=Value,colour=Color,shape=Stage),size=1.9,stroke=.55,fill='white',show.legend=TRUE) +
    geom_text(data=annotation,aes(x=center,y=y,label=Letters),family='Arial',fontface='bold',size=2.65,vjust=0,colour='#000000') +
    scale_color_manual(values=temp_colors,guide='none') +
    scale_fill_manual(values=temp_colors,guide='none') +
    scale_shape_manual(values=c(23,24,22),labels=c('Ancestral (T0)','Initial acclimation','Long-term selected (Tfinal)'),name=NULL,drop=FALSE) +
    scale_x_continuous(breaks=1:3,labels=c('23°C','28°C','32°C'),limits=c(.56,3.47),expand=expansion(0)) +
    labs(x=if(i==3) 'Measurement temperature' else NULL,y=if(j==1)ylabs[i] else NULL,
         title=if(i==1)fe_names[j] else ' ',tag=LETTERS[(i-1)*2+j]) +
    guides(shape=guide_legend(override.aes=list(colour=c('#676767','#000000','#000000'),size=2.4,stroke=.55,fill='white'),nrow=1)) +
    theme_classic(base_size=7,base_family='Arial') +
    theme(text=element_text(family='Arial',face='bold',colour='#000000'),
          axis.line=element_line(colour='#000000',linewidth=.60),axis.ticks=element_line(colour='#000000',linewidth=.50),
          axis.ticks.length=grid::unit(1.0,'mm'),axis.text=element_text(colour='#000000',size=7.5,face='bold'),
          axis.title.y=element_text(size=8.5,face='bold',colour='#000000',margin=margin(r=5)),axis.title.x=element_text(size=8.5,face='bold',colour='#000000',margin=margin(t=4)),
          plot.title=element_text(size=8.5,face='bold',colour='#000000',hjust=.5,margin=margin(b=4)),
          strip.background=element_blank(),strip.text=element_text(size=8.5,face='bold',colour='#000000',margin=margin(t=0,b=4)),
          panel.spacing.x=grid::unit(5,'mm'),legend.position='top',legend.text=element_text(size=7.5,colour='#000000',face='bold'),legend.title=element_text(size=7.5,colour='#000000',face='bold'),
          legend.key.width=grid::unit(4,'mm'),legend.spacing.x=grid::unit(3,'mm'),
          legend.margin=margin(0,0,3,0),plot.tag=element_text(size=9,face='bold',colour='#000000'),
          plot.tag.position='topleft',plot.margin=margin(t=3,r=5,b=3,l=3))
  # Freeze numeric formatting in each plot rather than retaining loop variable i.
  if(i==1) q <- q + scale_y_continuous(limits=c(0,limits[i]),breaks=breaks[[i]],expand=expansion(mult=c(0,.005)),labels=function(v) sprintf('%.1f',v))
  if(i>1) q <- q + scale_y_continuous(limits=c(0,limits[i]),breaks=breaks[[i]],expand=expansion(mult=c(0,.005)),labels=function(v) sprintf('%.2f',v))
  built <- ggplot_build(q); bx <- built$data[[2]]; pts <- built$data[[3]]
  stopifnot(nrow(pts)==nrow(dat), all(is.finite(pts$y)),all(abs(pts$y-dat$Value)<1e-12))
  stopifnot(nrow(built$data[[4]])==nrow(annotation), all(is.finite(built$data[[4]]$y)))
  rendered_letters[[length(rendered_letters)+1]] <- annotation[c('Trait','Treatment_Group','Fe','Letters','center','y')]
  for(g in unique(dat$Treatment_Group)) {
    gg <- dat[dat$Treatment_Group==g,];fe <- as.character(gg$Fe[1])
    row <- bx[abs(bx$x-gg$center[1])<1e-8,]
    stopifnot(nrow(row)==1, all(abs(gg$x-gg$center)<=.115))
    rendered_boxes[[length(rendered_boxes)+1]] <- data.frame(Treatment_Group=g,Trait=trait,Fe=fe,n=nrow(gg),
      Q1=row$lower,Median=row$middle,Q3=row$upper,Whisker_Lower=row$ymin,Whisker_Upper=row$ymax)
  }
  plots[[length(plots)+1]] <- q
}
write.csv(do.call(rbind,rendered_boxes),file.path(out,'source_data/Figure2_box_statistics.csv'),row.names=FALSE)
write.csv(do.call(rbind,rendered_letters),file.path(out,'source_data/Figure2_letter_positions.csv'),row.names=FALSE)
fig <- wrap_plots(G=guide_area(),A=plots[[1]],B=plots[[2]],C=plots[[3]],D=plots[[4]],E=plots[[5]],F=plots[[6]],
  design='GG\nAB\nCD\nEF',guides='collect',heights=c(.19,1.08,1,1.08)) +
  plot_annotation(theme=theme(plot.margin=margin(4,4,3,3)))
fig <- fig & theme(legend.position='top')
w <- 180/25.4; h <- 150/25.4
# Step 4: Save recreated figures and numerical drawing checks.
ggsave(file.path(out,'Figure2_revised.pdf'),fig,device=vector_pdf,width=w,height=h,family='Arial',bg='white')
ggsave(file.path(out,'Figure2_revised.svg'),fig,device=vector_svg,width=w,height=h,family='Arial',bg='white')
if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,'Figure2_revised.png'),fig,device=ragg::agg_png,width=w,height=h,units='in',dpi=300,bg='white')
if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(file.path(out,'Figure2_revised.tif'),fig,device=ragg::agg_tiff,width=w,height=h,units='in',dpi=600,compression='lzw',bg='white')
cat('Saved Figure 2 panels A-F with retained observations and one lowercase-letter family across ancestral, six initial-acclimation and six long-term-selected groups.\n')
