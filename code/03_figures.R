################################################################################
### Non-linearity and temporal variability are overlooked components of 
### global vertebrate population dynamics
### 04/2024

#### Part 3 - Figures
#### WARNING : Need to run file 02 before and have the outputs in the environments

################################################################################
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggpattern)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(tidyr)
library(treemapify)
library(sjstats) #cv function
library(forcats)
library(ggExtra)
library(moments)
library(lme4)
library(multcomp)
library(ggridges)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
#source("https://raw.githubusercontent.com/larmarange/JLutils/master/R/lm_right.R")

################################################################################
# Load the data -----------------------------------------------------------

final.nl <- read.csv("data/non_linear_models.csv")

################################################################################
# Functions ---------------------------------------------------------------

plot_classif_tot <- function(seed,data.long,data,col,xsp,ysp,lin=T){
  set.seed(seed)
  data<-data.long[data.long$ID==sample(data$ID,1),]
  sp<-stringr::str_replace(data$Binomial[1],"_"," ")
  abtype<-data.long$Units
  p<-ggplot(data=data,aes(x=year,y=logab))+
    geom_line(linewidth=1.5,color="grey50")
  if (lin){p<-p+geom_smooth(method="lm",se=FALSE,color=col,size=2)}
  else{p<-p+geom_smooth(method="lm",formula=y~poly(x,2),se=FALSE,color=col,size=2)}
  p<-p+theme_LPI()+ylab("log(abundance+1)")+
    annotate(geom="text", x=xsp, y=ysp, label=sp, size=5, family="Montserrat",fontface="italic")
}

theme_LPI <- function(){
  theme_classic() +
    theme(plot.title = element_text(size=14,family="Montserrat",face="bold",hjust=0.5),
          axis.text = element_text(size = 14,family="Montserrat"), 
          axis.title = element_text(size=14,family="Montserrat"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          # panel.grid.major.x = element_blank(),                                          
          # panel.grid.minor.x = element_blank(),
          # panel.grid.minor.y = element_blank(),
          # panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.7, 0.7, 0.7, 0.7), units = , "cm"))
}

drawWorld<-function(lats) {
  world_map<-map_data("world")
  g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
  g1<-g1+geom_polygon(data=world_map, aes(x=long, y=lat, group=group), colour="gray80", fill="gray80")
  g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
               panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
               axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
  return(g1)
}

plot_ex <- function(seed,data.long,data,lin=T){
  set.seed(seed)
  data<-data.long[data.long$ID==sample(data$ID,1),]
  p<-ggplot(data=data,aes(x=year,y=logab))
  if (lin){p<-p+geom_smooth(method="lm",se=FALSE,color="black",size=1.5)}
  else{p<-p+geom_smooth(method="lm",formula=y~poly(x,2),se=FALSE,color="black",size=1.5)}
  p<-p+theme_minimal()+theme(panel.grid.major = element_blank())
}

plot_prop <- function(p,col,line=F){
  p <- p + geom_bar(position = "fill") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("") +
    scale_fill_manual(values=col)
  if(line==T){p <- p+ geom_hline(yintercept = 0.448, linetype="dashed",color="grey30")}
  p <- p +  coord_flip() +
    theme_classic() +
    theme(text = element_text(family="Montserrat",face="bold"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank())
  p
}

plot_prop_traj <- function(p,col){
  p <- p + geom_bar(position = "fill") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("") +
    scale_fill_manual(values=col)
  #if(line==T){p <- p+ geom_hline(yintercept = sum(dframe_taxa$LoNL=="non linear")/6437, linetype="dashed",color="grey30")}
  p <- p +  coord_flip() +
    theme_classic() +
    theme(text = element_text(family="Montserrat",face="bold"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank())
  p
}

plot_prop2 <- function(p,col,line=F){
  p <- p + geom_bar(position = "fill") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("") +
    scale_fill_manual(values=col)
  if(line==T){p <- p+ geom_hline(yintercept = sum(dframe_taxa$LoNL=="non linear")/6437, linetype="dashed",color="grey30")}
  p <- p +  
    theme_classic() +
    theme(text = element_text(family="Montserrat",face="bold"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank())
  p
}

plot_number <- function(df,ctgr,col){
  p <- ggplot(df, aes(x = fct_rev(fct_infreq(ctgr)), fill = ctgr)) +
     geom_bar() +
     geom_text(aes(label = after_stat(count)), hjust=1,stat = "count",family="Montserrat",fontface="bold",color="black") +
     xlab("") + ylab("Number of population") +
     scale_fill_manual(values=col)+
     coord_flip()+ theme_classic() +
     theme(text = element_text(family="Montserrat",face="bold"),
           legend.position = "none",
           axis.ticks.y = element_blank(),
           axis.line.y = element_blank())
  p
}

################################################################################
# Colors ------------------------------------------------------------------
cols <- c("#8F0606", "#D1AB41", "#426962")
cols22 <- c("#8F0606", "#426962", "#D1AB41")
cols2 <- c("#6E2E2E", "#AD8080", "#4D806A", "#739989", "#8F7232", "#E0BD72")
colstraj <- c("#A84646","#C4993D","#4D806A","#739989","#E0BD72","#AD8080")
colprop<-c("#4c8cac", "#cf994c")
cols3 <- c("#156085", "#FF9721","#156085", "#FF9721","#156085", "#FF9721")
cols4 <- c("#A84646","#6E2E2E","#AD8080",
           "#4D806A","#38574A", "#739989",
           "#C4993D","#E0BD72","#8F7232")
cols5 <- c("#6E2E2E","#A84646","#38574A",
           "#4D806A","#C4993D","#8F7232",
           "#AD8080", "#739989","#E0BD72")
colssyst <- c("#4D806A","#273045","#8F7232")
colssyst <- c("#0b775e", "#273046", "#a2a475")
colstaxa <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1")
colsregion <- c("#AD8080","#72B5DB","#4D806A","#8F7232","#156085","#FF9721","#AA9AB8","#FFD16E")
colsrealms=c("#AD8080",'#94b6d2','#704404','#34405e',
             '#d8b25c','#f7b615','#968c8c','#dd8047',
             '#ebddc3','#7ba79d','#a5ab81','#38784d',
             '#775f55')
colstaxawhole <- c("#4D806A","#2FA1A1","#8F7232","#FF9721","#156085","#AD8080",
                   "#4D806A","#2FA1A1","#8F7232","#FF9721","#156085","#AD8080")
colstaxa2 <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1",
               "#8F7232","#156085","#AD8080","#FF9721","#2FA1A1")
colstaxa3 <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1",
               "#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1")
colstaxa3 <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1",
               "darkseagreen","#C9AF72","skyblue3","rosybrown2","#FFBC70","#88D4D4")
colstaxa4 <- c("darkseagreen","#C9AF72","skyblue3","rosybrown2","#FFBC70","#88D4D4",
               "#C9AF72","skyblue3","rosybrown2","#FFBC70","#88D4D4")

################################################################################
############################### MAIN FIGURES ###################################
################################################################################

################################################################################
# Figure 2 | All classifications ------------------------------------------

# A) Linear decrease
lin_dec <- final.nl  %>% filter(shape_class=="decrease_constant" & lin_class=="negative")
lin_dec.long <- lpd.long %>% dplyr::filter(ID %in% lin_dec$ID)
(A<-plot_classif_tot(8,lin_dec.long,lin_dec,col="#A84646",xsp=1995,ysp=0.55)+
    annotate(geom="text", x=2011, y=0.78, label="N=1170\n18.2%", family="Montserrat",fontface="bold",size=5,col="#A84646")+
    ggtitle("Linear decrease"))

# B) No trend (linear)
no_trend <- final.nl  %>% filter(shape_class=="stable_constant" & lin_class=="no trend")
no_trend.long <- lpd.long %>% dplyr::filter(ID%in% no_trend$ID)
(B<-plot_classif_tot(33,no_trend.long,no_trend,col="#C4993D",xsp=1985,ysp=1)+
    annotate(geom="text", x=2004.5, y=1.68, label="N=1337\n20.8%", family="Montserrat",fontface="bold",size=5,col="#C4993D")+
    ggtitle("No trend"))

# C) Linear increase
lin_inc <- final.nl  %>% filter(shape_class=="increase_constant" & lin_class=="positive")
lin_inc.long <- lpd.long %>% dplyr::filter(ID%in% lin_inc$ID)
(C<-plot_classif_tot(13,lin_inc.long,lin_inc,col="#4D806A",xsp=2005,ysp=4.45)+
    annotate(geom="text", x=1980, y=4.93, label="N=1043\n16.2%", family="Montserrat",fontface="bold",size=5,col="#4D806A")+
    ggtitle("Linear increase"))

# D) Accelerated decrease
accel_dec <- final.nl %>% filter(shape_class=="decrease_accelerated" & lin_class=="negative")
accel_dec.long <- lpd.long %>% dplyr::filter(ID%in% accel_dec$ID)
(D<-plot_classif_tot(42,accel_dec.long,accel_dec,col="#6E2E2E",xsp=1982,ysp=2.2,lin=F)+
    annotate(geom="text", x=2010, y=2.76, label="N=215\n3.3%", family="Montserrat",fontface="bold",size=5,col="#6E2E2E")+
    ggtitle("Decrease accelerated"))

# E) Concave
concave <- final.nl  %>% filter(shape_class=="stable_concave" & lin_class=="no trend")
concave.long <- lpd.long %>% dplyr::filter(ID%in% concave$ID)
(E<-plot_classif_tot(4,concave.long,concave,col="#E0BD72",xsp=1990,ysp=4.5,lin=F)+
    annotate(geom="text", x=1993, y=8.5, label="N=757\n11.8%", family="Montserrat",fontface="bold",size=5,col="#E0BD72")+
    ggtitle("Concave"))

# F) Decelerated increase
decel_inc <- final.nl  %>% filter(shape_class=="increase_decelerated" & lin_class=="positive")
decel_inc.long <- lpd.long %>% dplyr::filter(ID%in% decel_inc$ID)
(Fi<-plot_classif_tot(6,decel_inc.long,decel_inc,col="#739989",xsp=1985,ysp=6.7,lin=F)+
    annotate(geom="text", x=1955, y=8.75, label="N=276\n4.3%", family="Montserrat",fontface="bold",size=5,col="#739989")+
    ggtitle("Increase decelerated"))

# G) Decelerated decrease
decel_dec <- final.nl %>% filter(shape_class=="decrease_decelerated" & lin_class=="negative")
decel_dec.long <- lpd.long %>% dplyr::filter(ID%in% decel_dec$ID)
(G<-plot_classif_tot(57,decel_dec.long,decel_dec,col="#AD8080",xsp=1965,ysp=6.4,lin=F)+
    annotate(geom="text", x=2006, y=9, label="N=554\n8.6%", family="Montserrat",fontface="bold",size=5,col="#AD8080")+
    ggtitle("Decrease decelerated"))

# H) Convex
convex <- final.nl  %>% filter(shape_class=="stable_convex" & lin_class=="no trend")
convex.long <- lpd.long %>% dplyr::filter(ID%in% convex$ID)
(H<-plot_classif_tot(15,convex.long,convex,col="#8F7232",xsp=1972,ysp=5.7,lin=F)+
    annotate(geom="text", x=1997.5, y=5.8, label="N=621\n9.6%", family="Montserrat",fontface="bold",size=5,col="#8F7232")+
    ggtitle("Convex"))

# I) Accelerated increase
accel_inc <- final.nl  %>% filter(shape_class=="increase_accelerated" & lin_class=="positive")
accel_inc.long <- lpd.long %>% dplyr::filter(ID%in% accel_inc$ID)
(I<-plot_classif_tot(61,accel_inc.long,accel_inc,col="#38574A",xsp=1989,ysp=0,lin=F)+
    annotate(geom="text", x=1972, y=1.15, label="N=464\n7.2%", family="Montserrat",fontface="bold",size=5,col="#38574A")+
    ggtitle("Increase accelerated"))

(all<-cowplot::plot_grid(A,B,C,D,E,Fi,G,H,I,labels=c("A","B","C","D","E","F","G","H","I")))
ggsave("outputs/all_classif.png",all,type='cairo',width=15,height=12,dpi=200)

################################################################################
# Figure 3 | Cross-classification -----------------------------------------
final.nl$shape_class <- as.factor(final.nl$shape_class)

double_classif <- final.nl %>% 
  mutate(shape_class=fct_relevel(shape_class,
                                 "decrease_constant","decrease_accelerated","decrease_decelerated",
                                 "increase_constant","increase_accelerated","increase_decelerated",  
                                 "stable_constant","stable_concave","stable_convex"))

(cross_classif <- ggplot(double_classif, aes(x = lin_class, fill = shape_class)) +
    geom_bar(position = "fill",width = 0.25) +
    geom_text(aes(label = ifelse(after_stat(count)>10,after_stat(count),"")), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("Proportion") +
    scale_fill_manual(values=cols4)+ theme_void() + labs(fill = "Linearity") + 
    theme(text = element_text(family="Montserrat",face="bold"),legend.position = "none"))

ggsave("outputs/cross_classif.png",cross_classif,width = 7, height = 6)


################################################################################
# Figure 4 | Non-linearity across habitat types, regions ------------------

plot_prop <- function(x){
  plot(x)+xlab("")+ylab("")+xlim(c(0.1,0.9))+geom_vline(xintercept = 0.5,linetype="dashed")+
    theme_bw() + theme(text = element_text(family="Montserrat",face="bold"),
                       axis.ticks.y = element_blank())
}

plot_prop_var <- function(x){
  plot(x)+
    xlab("")+ylab("")+xlim(c(0,NA))+
    #geom_vline(xintercept = mean(d$MSE),linetype="dashed")+
    theme_bw() + theme(text = element_text(family="Montserrat",face="bold"),
                       axis.ticks.y = element_blank())
}

# Habitat
plot_hab <- plot_prop(EMM_hab)+ggtitle("A. Habitat type")
#geom_text(data = hab.let_df, aes(label = letters, y=System, x = 0.9 ), size=4, family="Montserrat",fontface="italic",colour="black", size=5)

# Region
plot_reg <- plot_prop(EMM_reg)+ggtitle("B. Region")+  scale_y_discrete(labels = scales::label_wrap(20)) +
  geom_text(data = reg.let_df, aes(label = letters, y=Region, x = 0.9 ), size=4, family="Montserrat",fontface="italic",colour="black", size=5)

# Taxonomic group
plot_taxa <- plot_prop(EMM_taxa)+ggtitle("C. Taxonomic group")+
  geom_text(data = taxa.let_df, aes(label = letters, y=taxa, x = 0.9 ), size=4, family="Montserrat",fontface="italic",colour="black", size=5)

# RLC
plot_iucn <- plot_prop(EMM_iucn)+xlab("Proportion of non-linearity (adjusted after GLMM fit)")+ggtitle("D. IUCN Red List Category")
#geom_text(data = rlc.let_df, aes(label = letters, y=redlistCategory, x = 0.9 ), size=4, family="Montserrat",fontface="italic",colour="black", size=5)

# All
all_lin <- cowplot::plot_grid(plot_hab,plot_reg,plot_taxa,plot_iucn,ncol=1,rel_heights = c(4,9,7,7),align="hv")
all_lin

ggsave("outputs/all_lin.png",all_lin,height = 10,width=8)

################################################################################
# Figure 5 | Temporal variability across habitat types, regions -----------

(plot_hab <- plot_prop_var(EMM_full_hab)+ggtitle("A.Habitat type")+
    geom_text(data = hab.let_full_df, aes(label = letters, y=System, x = 0.3 ), 
              size=4, family="Montserrat",fontface="italic",colour="black", size=5))

(plot_reg <- plot_prop_var(EMM_full_reg)+ggtitle("B. Region")+
    scale_y_discrete(labels = scales::label_wrap(20)) +
    geom_text(data = reg.let_full_df, aes(label = letters, y=Region, x = 0.3 ), 
              size=4, family="Montserrat",fontface="italic",colour="black", size=5))

(plot_taxa <- plot_prop_var(EMM_full_taxa)+ggtitle("C. Taxonomic group")+
    geom_text(data = taxa.let_full_df, aes(label = letters, y=taxa, x = 0.3 ), 
              size=4, family="Montserrat",fontface="italic",colour="black", size=5))

(plot_iucn <- plot_prop_var(EMM_full_iucn)+xlab("MSE (adjusted after GLMM fit)")+
    ylab("")+ggtitle("D. IUCN Red List Category")+
    geom_text(data = rlc.let_df, aes(label = letters, y=redlistCategory, x = 5e+07 ), 
              size=4, family="Montserrat",fontface="italic",colour="black", size=5))

(all_var <- cowplot::plot_grid(plot_hab,plot_reg,plot_taxa,plot_iucn,ncol=1,
                               rel_heights = c(4,9,7,7),align="hv"))

ggsave("outputs/all_var.png",all_var,height = 10,width=8)


################################################################################
########################## SUPPLEMENTARY FIGURES ###############################
################################################################################

################################################################################
# SM1 | Temporal, geographical, and taxonomic extent ----------------------

### FS1.1 | DURATION 
(duration.all <- ggplot(final.nl, aes(points)) +
   geom_histogram(binwidth = 1, alpha = 0.2, position = "identity") +
   geom_line(stat = "density", aes(y = ..count..), size = 1.5) +
   theme_classic() +
   scale_y_continuous(limits = c(0, 450)) +
   geom_vline(xintercept = mean(final.nl$points), size = 1) +
   scale_x_continuous(breaks = c(20, 25, 30, 35, 40, 45, 50, 55,60,65,70)) +
   theme(axis.line.x = element_line(color = "black"),
         axis.line.y = element_line(color = "black"),
         text = element_text(family="Montserrat",size=15,face="bold")) +
   labs(x = "Number of monitoring years", y = "Number of time series") +
   geom_vline(xintercept = 19.5, colour = "grey30", size = 1, linetype = "dashed") +
   guides(colour = F, fill = F))
ggsave("outputs/temporal_extent.png",duration.all,width=15,height = 10)

### FS1.2 | DISTRIBUTIONS ACROSS BIOGEOGRAPHIC AND TAXONOMIC GROUPS
barhab<-plot_number(final.nl,final.nl$System,colssyst)
bartax<-plot_number(final.nl %>% filter(taxa!="Invertebrates"),(final.nl %>% filter(taxa!="Invertebrates"))$taxa,colstaxa)
barreg<-plot_number(final.nl,final.nl$Region,colsregion)
barrealm<-plot_number(final.nl,final.nl$full_realm,colsrealms)
bariucn<-plot_number(data.rlc,data.rlc$redlistCategory,c("#ff000e","#ffa15b","#ffef00", "#d5db00","#4bb33c",'grey'))

co11<-cowplot::plot_grid(NULL,barhab,ncol=2,rel_widths = c(0.2,1))
co21<-cowplot::plot_grid(NULL,bartax,ncol=2,rel_widths = c(0.15,1))
co1<-cowplot::plot_grid(co11,barreg,labels=c("A","C"),ncol=1,rel_heights = c(0.3,0.8))
co2<-cowplot::plot_grid(co21,bariucn,labels=c("B","D"),ncol=1,rel_heights = c(0.6,0.7))
co3<-cowplot::plot_grid(co1,co2,ncol=2)
(co4<-cowplot::plot_grid(co3,barrealm,labels=c("","E"),ncol=1))
ggsave("outputs/extents.png",co4,height = 12,width=15)

### FS1.3 | MAPS
(world.map.hab <- drawWorld() +
   geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=System),size=1,alpha = 0.7)+ 
   scale_colour_manual(values=colssyst) +
   scale_y_continuous(limits = c(-80, 80)) +
   #theme(legend.position = "bottom") +
   guides(color = guide_legend(override.aes = list(size = 4,alpha=1))))

(world.map.tax <- drawWorld() +
    geom_point(data = final.nl %>% filter(taxa!="Invertebrates"), aes(x = Longitude, y = Latitude, color=taxa),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=colstaxa) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.reg <- drawWorld() +
    geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=Region),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=colsregion) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.realm <- drawWorld() +
    geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=full_realm),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=colsrealms) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.iucn <- drawWorld() +
    geom_point(data = data.rlc, aes(x = Longitude, y = Latitude, color=redlistCategory),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=c("#ff000e","#ffa15b","#ffef00", "#d5db00","#4bb33c",'grey')) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.nl <- drawWorld() +
    geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=LoNL),alpha = 0.7)+ 
    scale_colour_manual(values=colprop) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size = guide_legend(title = "Number of years monitored"),
           color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

# ggsave(world.map.nl, filename = "world.map.nl.pdf", device = "pdf", dpi = 300,
#        height = 15, width = 15)

(maxiplot<-cowplot::plot_grid(world.map.hab,
                              world.map.tax,
                              world.map.reg,
                              world.map.realm,
                              world.map.iucn,
                              world.map.nl,
                              ncol=2, labels=c("A","B","C","D","E","F")))

ggsave("outputs/maps.png",maxiplot,height = 12,width=20)

### TS1.1 | CROSS-DISTRIBUTION HABITAT TYPES AND REGIONS
table(final.nl$Region,final.nl$System)

################################################################################
# SM3 | Impact of the duration, number of years samples, and start --------

test_firstX <- final.nl %>% 
  group_by(first_X_value) %>% 
  summarise(prop=sum(LoNL=="non linear")/n(),
            n=n())

test_points <- final.nl %>% 
  group_by(points) %>% 
  summarise(prop=sum(LoNL=="non linear")/n(),
            n=n())

p1supp <- ggplot(test_firstX,aes(x=first_X_value,y=prop))+geom_line()+
  geom_smooth(method="lm",color="black")+ theme_classic()+
  xlab("Starting year of time series") + ylab("Proportion of non-linearity") +
  theme(text = element_text(family="Montserrat",face="bold"),legend.position = "bottom")

p2supp <- ggplot(test_points,aes(x=points,y=prop))+geom_line()+
  geom_smooth(method="lm",color="black")+ theme_classic()+
  xlab("Number of points") + ylab("Proportion of non-linearity") +
  theme(text = element_text(family="Montserrat",face="bold"),legend.position = "bottom")

psupp <- cowplot::plot_grid(p1supp,p2supp,ncol=2,labels=c("A","B"))

p3supp <- ggplot(final.nl, aes(x = first_X_value, y = points5cl, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1, bandwidth=2) +
  scale_fill_viridis(name = "first_X_value", option = "C") + theme_classic() +
  theme(legend.position = "none", text = element_text(family="Monserrat",face="bold",size=12))+
  xlab("Starting year of the time series")+ylab("Number of points within the time series")


supp<-cowplot::plot_grid(psupp,p3supp,ncol=1,labels=c("","C"),rel_heights = c(1,1.5))
supp
ggsave("outputs/side_analysis.png",supp,height = 10,width=10)

################################################################################
# SM4 | Investigating the sensitivity of the results using differe --------

plot_hab_D <- plot_prop_var(EMM_hab_D)+ggtitle("A.Habitat type")+
  geom_text(data = hab.let_D_df, aes(label = letters, y=System, x = 0.55 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)

plot_reg_D <- plot_prop_var(EMM_reg_D)+ggtitle("C. Region")+
  scale_y_discrete(labels = scales::label_wrap(20)) +
  geom_text(data = reg.let_D_df, aes(label = letters, y=Region, x = 0.55 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)

plot_taxa_D <- plot_prop_var(EMM_taxa_D)+ggtitle("E. Taxonomic group")+
  geom_text(data = taxa.let_D_df, aes(label = letters, y=taxa, x = 0.55 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)

plot_iucn_D <- plot_prop_var(EMM_iucn_D)+xlab("D index (adjusted after GLMM fit)")+
  ylab("")+ggtitle("G. IUCN Red List Category")

all_var_D <- cowplot::plot_grid(plot_hab_D,plot_reg_D,plot_taxa_D,plot_iucn_D,ncol=1,
                                rel_heights = c(4,9,7,7),align="hv")

plot_hab_CV <- plot_prop_var(EMM_hab_CV)+ggtitle("B.Habitat type")+
  geom_text(data = hab.let_CV_df, aes(label = letters, y=System, x = 0.8 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)+
  theme(axis.text.y = element_blank())

plot_reg_CV <- plot_prop_var(EMM_reg_CV)+ggtitle("D. Region")+
  scale_y_discrete(labels = scales::label_wrap(20)) +
  geom_text(data = reg.let_CV_df, aes(label = letters, y=Region, x = 1.2 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)+
  theme(axis.text.y = element_blank())

plot_taxa_CV <- plot_prop_var(EMM_taxa_CV)+ggtitle("F. Taxonomic group")+
  geom_text(data = taxa.let_CV_df, aes(label = letters, y=taxa, x = 1 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)+
  theme(axis.text.y = element_blank())

plot_iucn_CV <- plot_prop_var(EMM_iucn_CV)+xlab("CV (adjusted after GLMM fit)")+
  ylab("")+ggtitle("H. IUCN Red List Category")+
  theme(axis.text.y = element_blank())

all_var_CV <- cowplot::plot_grid(plot_hab_CV,plot_reg_CV,plot_taxa_CV,plot_iucn_CV,ncol=1,
                                 rel_heights = c(4,9,7,7),align="hv")

(all_var_supp<-cowplot::plot_grid(all_var_D,all_var_CV,ncol=2,rel_widths = c(2,1.6)))
ggsave("outputs/all_var_supp.png",all_var_supp,height = 10,width=10)


################################################################################
# SM5 | Detailed repartition of non-linear trajectories among biog --------

### Donuts Taxonomic groups
taxa_df <- as.data.frame(table(final.nl$taxa,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var1=fct_relevel(Var1,c("Amphibians","Invertebrates","Sharks_Rays","Birds","Reptiles","Fish","Mammals")),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(taxa_donuts <- ggplot(data=taxa_df %>% filter(Var1!="Invertebrates"), 
                       aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,ncol=2) +theme_void()+ 
    theme(legend.position = "none", 
          # legend.box = "vertical",
          text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    guides(fill=guide_legend(title="Type of trajectory",nrow=1)) +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

### Donuts Habitat types
hab_df <- as.data.frame(table(final.nl$System,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(hab_donuts <- ggplot(data=hab_df, aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,nrow=1) +theme_void()+ 
    theme(legend.position = "none", text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

### Donuts Regions
reg_df <- as.data.frame(table(final.nl$Region,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(regions_donuts <- ggplot(data=reg_df, aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,nrow=2) +theme_void()+ 
    theme(legend.position = "bottom", legend.box = "vertical",
          text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    
    guides(fill=guide_legend(title="",nrow=1),color="none") +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

### Donuts RLC
rlc_df <- as.data.frame(table(final.nl$redlistCategory,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var1=fct_relevel(Var1,c("Extinct in the Wild","Critically Endangered","Endangered",
                                 "Vulnerable","Near Threatened","Least Concern","Data Deficient")),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(rlc_donuts <- ggplot(data=rlc_df %>% filter(Var1!="Extinct in the Wild"), 
                      aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,ncol=2) +theme_void()+ 
    theme(legend.position = "none", legend.box = "vertical",
          text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    guides(fill=guide_legend(title="",nrow=1),color="none") +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

(codoreg <- cowplot::plot_grid(hab_donuts,regions_donuts,labels=c("A","B"),ncol=1,rel_heights = c(0.75,2)))
(codotax <- cowplot::plot_grid(taxa_donuts,rlc_donuts,labels=c("A","B"),ncol=2))
ggsave("outputs/donuts_reg.png",codoreg,height=10,width=10)
ggsave("outputs/donuts_tax.png",codotax,height=10,width=12)


################################################################################
# SM6 | Temporal variability across trajectory types

SM6A <- plot_prop_var(EMM_full)+
  geom_text(data = traj.let_df, aes(label = letters, y=shape_class, x = 0.3 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)+
  xlab("MSE (adjusted after GLMM fit)")+ggtitle("A. MSE")

SM6B <- plot_prop_var(EMM_full_D)+
  geom_text(data = traj.let_D_df, aes(label = letters, y=shape_class, x = 0.5 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)+
  xlab("D index (adjusted after GLMM fit)")+ggtitle("B. D index")

SM6C <- plot_prop_var(EMM_full_CV)+
  geom_text(data = traj.let_CV_df, aes(label = letters, y=shape_class, x = 0.4 ), 
            size=4, family="Montserrat",fontface="italic",colour="black", size=5)+
  xlab("CV (adjusted after GLMM fit)")+ggtitle("C. CV")

(trajvar <- cowplot::plot_grid(SM6A,SM6B,SM6C,ncol=1))
ggsave("outputs/trajvar.png",trajvar,height=10,width=8)
