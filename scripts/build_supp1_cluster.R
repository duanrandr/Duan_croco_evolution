#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
# Step 1: Locate repository inputs and create a separate output directory.
code_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE))))
repo <- normalizePath(file.path(code_dir, '..', '..'))
out <- file.path(repo, 'results', 'reproduced', 'multivariate')
dir.create(file.path(out, 'source_data'), recursive=TRUE, showWarnings=FALSE)
src <- file.path(repo, 'data', 'figure_source_data', 'multivariate')
# Step 2: Load export devices and read the supplied figure inputs.
source(file.path(code_dir, 'vector_export.R'))
src <- file.path(out, 'source_data')
z <- read.csv(file.path(src, 'Standardized_matrix.csv'), stringsAsFactors=FALSE, check.names=FALSE)
link <- read.csv(file.path(src, 'Average_linkage.csv'), stringsAsFactors=FALSE)
leaf <- read.csv(file.path(src, 'Leaf_order_bottom_to_top.csv'), stringsAsFactors=FALSE)
traits <- c('Growth_rates','Nfix_PON','Cfix_POC')
metadata <- read.csv(file.path(src,'Treatment_mean_matrix.csv'),stringsAsFactors=FALSE)
target <- metadata$Treatment_Group[as.logical(metadata$Highlight)]
n <- nrow(z)
stopifnot(n==25, nrow(link)==24, nrow(leaf)==25,
          all(link$Node_ID==n:(2*n-2)),
          setequal(leaf$Treatment_Group,z$Treatment_Group),
          all(z$Treatment_Group[leaf$Input_Index+1]==leaf$Treatment_Group),
          all(is.finite(as.matrix(z[traits]))))

label_only <- FALSE  # Step 3: Rebuild and validate the complete supplied tree.
if(label_only) {
  # For wording-only updates, reuse the exact validated drawing coordinates.
  # No clustering, standardization, or source-data write is performed.
  edges <- read.csv(file.path(src,'SuppFig1_rendered_tree_segments.csv'),stringsAsFactors=FALSE)
  validation <- read.csv(file.path(src,'SuppFig1_linkage_validation.csv'),stringsAsFactors=FALSE)
  leaf <- read.csv(file.path(src,'SuppFig1_rendered_leaf_labels.csv'),stringsAsFactors=FALSE)
  stopifnot(nrow(edges)==3*(n-1),nrow(leaf)==n,
            setequal(leaf$Treatment_Group,z$Treatment_Group))
} else {
# Independently validate every supplied merge against R's average-linkage tree.
# Tree drawing below uses the supplied linkage and original leaf order exactly.
hc <- hclust(dist(z[traits],method='euclidean'),method='average')
r_members <- vector('list',n-1)
for(i in seq_len(n-1)) {
  get_child <- function(k) if(k<0) -k else r_members[[k]]
  r_members[[i]] <- sort(c(get_child(hc$merge[i,1]),get_child(hc$merge[i,2])))
}
signature <- function(v) paste(sort(v),collapse=',')
r_keys <- vapply(r_members,signature,character(1))
members <- vector('list',2*n-1)
node_x <- rep(NA_real_,2*n-1); node_y <- rep(NA_real_,2*n-1)
pure <- rep(FALSE,2*n-1)
for(i in seq_len(n)) {
  members[[i]] <- i
  node_x[i] <- 0
  node_y[i] <- leaf$Position[match(i-1,leaf$Input_Index)]
  pure[i] <- z$Treatment_Group[i] %in% target
}
edges <- list(); validation <- list()
for(i in seq_len(nrow(link))) {
  p <- link$Node_ID[i]+1; a <- link$Child_1[i]+1; b <- link$Child_2[i]+1
  stopifnot(a<p,b<p)
  members[[p]] <- sort(c(members[[a]],members[[b]]))
  node_x[p] <- link$Height[i]
  node_y[p] <- mean(c(node_y[a],node_y[b]))
  pure[p] <- pure[a] && pure[b]
  r_id <- match(signature(members[[p]]),r_keys)
  stopifnot(!is.na(r_id),length(members[[p]])==link$Count[i],
            abs(hc$height[r_id]-link$Height[i])<1e-12)
  validation[[i]] <- data.frame(Node_ID=link$Node_ID[i],Count=link$Count[i],
    Supplied_Height=link$Height[i],R_Height=hc$height[r_id],
    Absolute_Difference=abs(hc$height[r_id]-link$Height[i]),
    Members=paste(z$Treatment_Group[members[[p]]],collapse=';'))
  # Highlight only merges contained entirely within the three target treatments.
  edges[[i]] <- data.frame(Node_ID=link$Node_ID[i],Highlighted=pure[p],
    x=c(node_x[a],node_x[b],node_x[p]),
    xend=rep(node_x[p],3),
    y=c(node_y[a],node_y[b],node_y[a]),
    yend=c(node_y[a],node_y[b],node_y[b]))
}
edges <- do.call(rbind,edges)
validation <- do.call(rbind,validation)
stopifnot(all(target %in% z$Treatment_Group))
write.csv(validation,file.path(src,'SuppFig1_linkage_validation.csv'),row.names=FALSE)
write.csv(edges,file.path(src,'SuppFig1_rendered_tree_segments.csv'),row.names=FALSE)
}

leaf$Label <- metadata$Dendrogram_Label[match(leaf$Treatment_Group,metadata$Treatment_Group)]
leaf$Label <- sub('^Initial acclimation', 'IA', leaf$Label)
leaf$Highlighted <- leaf$Treatment_Group %in% target
leaf$Face <- 'bold'
leaf$Colour <- '#000000'
if(!label_only) write.csv(leaf,file.path(src,'SuppFig1_rendered_leaf_labels.csv'),row.names=FALSE)
if(label_only) write.csv(leaf,file.path(out,'qa/IA_label_only/final_rendered_leaf_labels.csv'),row.names=FALSE)

# Preserve the original brown emphasis only on pure clades actually recovered.
# The membership and tree topology are never forced to keep the old highlight.
if(label_only) {
  highlight_positions <- sort(leaf$Position[leaf$Highlighted])
} else {
  highlight_members <- unique(unlist(members[which(pure & lengths(members)>=2)]))
  highlight_positions <- sort(leaf$Position[(leaf$Input_Index+1) %in% highlight_members])
}
bands <- data.frame(ymin=numeric(),ymax=numeric())
if(length(highlight_positions)>0) {
  runs <- split(highlight_positions,cumsum(c(TRUE,diff(highlight_positions)!=1)))
  bands <- do.call(rbind,lapply(runs,function(p)data.frame(ymin=min(p)-.48,ymax=max(p)+.48)))
}
xmax <- max(3.5,ceiling(max(link$Height)*2)/2)

fig <- ggplot() +
  geom_rect(data=bands,aes(ymin=ymin,ymax=ymax),xmin=-1.77,xmax=xmax,
    fill='#f3eae6',colour=NA) +
  geom_segment(data=edges[!edges$Highlighted,],
    aes(x=x,y=y,xend=xend,yend=yend),colour='#808080',linewidth=.43,lineend='square') +
  geom_segment(data=edges[edges$Highlighted,],
    aes(x=x,y=y,xend=xend,yend=yend),colour='#4a1f1f',linewidth=.85,lineend='square') +
  geom_text(data=leaf,aes(x=-.09,y=Position,label=Label,fontface=Face,colour=Colour),
    hjust=1,family='Arial',size=7.5/ggplot2::.pt,show.legend=FALSE) +
  scale_colour_identity() +
  scale_x_continuous(breaks=seq(0,xmax,.5),labels=function(x)sprintf('%.1f',x),expand=expansion(0)) +
  scale_y_continuous(breaks=NULL,expand=expansion(0)) +
  coord_cartesian(xlim=c(0,xmax),ylim=c(.40,n+.65),clip='off') +
  labs(x='Euclidean distance (standardized traits)',y=NULL) +
  theme_classic(base_size=8,base_family='Arial') +
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
    axis.line.x=element_line(linewidth=.60,colour='#000000'),
    axis.ticks.x=element_line(linewidth=.50,colour='#000000'),
    axis.ticks.length=grid::unit(1.25,'mm'),
    axis.text.x=element_text(size=7.5,face='bold',colour='#000000'),
    axis.title.x=element_text(size=8.5,face='bold',colour='#000000',margin=margin(t=5)),
    plot.margin=margin(3,3,3,62,unit='mm'))

built <- ggplot_build(fig)
stopifnot(nrow(built$data[[2]])+nrow(built$data[[3]])==3*(n-1),
          nrow(built$data[[4]])==n,all(is.finite(built$data[[4]]$y)))
# Step 4: Save recreated figures and numerical drawing checks.
stem <- file.path(out,'SuppFig1_hierarchical_clustering')
w <- 180/25.4; h <- 145/25.4
ggsave(paste0(stem,'.pdf'),fig,device=vector_pdf,width=w,height=h,family='Arial',bg='white')
ggsave(paste0(stem,'.svg'),fig,device=vector_svg,width=w,height=h,family='Arial',bg='white')
if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(paste0(stem,'.png'),fig,device=ragg::agg_png,width=w,height=h,dpi=300,bg='white')
if(Sys.getenv('CROCO_VECTOR_ONLY') != '1') ggsave(paste0(stem,'.tif'),fig,device=ragg::agg_tiff,width=w,height=h,dpi=600,compression='lzw',bg='white')
if(label_only) {
  cat('Saved label-only SuppFig1 from existing tree coordinates; no analysis data were recalculated or written.\n')
} else {
  cat(sprintf('Saved SuppFig1; all %d merge clades match R average linkage; maximum height difference %.3g.\n',n-1,max(validation$Absolute_Difference)))
}
