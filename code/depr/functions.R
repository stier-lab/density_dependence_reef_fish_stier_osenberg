#relevant functions for individual study analysis and bootstrapping meta analysis


#####
#Perty figs
#####

#wes anderson color palette 
pal <- wes_palette("Zissou1", 100, type = "continuous")

gc <- guide_colorbar(
  frame.colour = "black",
  barheight = 8,
  frame.linewidth = 2,
  ticks.colour = "black",
  ticks.linewidth = 2
)

#bars to background of plot
# geom_stripes(
#   mapping = NULL,
#   data = NULL,
#   stat = "identity",
#   position = "identity",
#   ...,
#   show.legend = NA,
#   inherit.aes = TRUE
# )


#write down inverse hyperbolic sin transformation
# https://robjhyndman.com/hyndsight/transformations/
# https://stackoverflow.com/questions/14504869/histogram-with-negative-logarithmic-scale-in-r
# https://stats.stackexchange.com/questions/1444/how-should-i-transform-non-negative-data-including-zeros

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}


limits <- 100
step <- 0.005
demo <- data.frame(x=seq(from=-1*limits,to=limits,by=step))


ggplot(demo,aes(x,x))+geom_point(size=2)+
  scale_y_continuous(trans = 'asinh',breaks=c(-100,-50,-10,-1,0,1,10,50,100))+
  theme_bw()


#####
#Try Catch - Loop without stopping
#####

tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler), warning = W)
}


####
# Model plotting functions 
####

bh_f = function(alpha,beta,t_days,x) {
  recruits = (x*(exp(-alpha*t_days)))/(1+(beta*x*(1-exp(-alpha*t_days))/alpha))
  return(recruits)
}

bh_prob = function(settlers) {
  (exp(-alpha*t_days))/(1+(beta*settlers*(1-exp(-alpha*t_days))/alpha))
}



###################
### Bootstrapping Functions
###################

# r1boot <- function() {
#   bdat <- dat[sample(nrow(dat),size=nrow(dat),replace=TRUE),]  ## bootstrap sample
#   rb <- rma(yi=betamle2, vi=beta.semle2, data=dat, method="REML",test="knha")      ## re-run model
#   as.matrix(coef(rb))                                          ## extract/reformat results
# }

#prefer this one http://www.metafor-project.org/doku.php/tips:bootstrapping_with_ma non parametirc bootstrap to see if it differes or how


#non-parametric bootstrapping function from metafor website no random effects
boot.func <- function(dat, indices) {
  
  sel <- dat[indices,]
  res <- try(rma(yi=betanls2, vi=betanlsvar, data=sel, subset=indices,method="REML", test="knha"), silent=TRUE)
  if (is.element("try-error", class(res))) {
    NA
  } else {
    c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
}

#non-parametric bootstrapping function from metafor website adapted to use for mixed model with study/substudy
boot.func_2 <- function(dat, indices) {
  
  res <- try(rma.mv(yi=betanls2, V=betanlsvar, data=dat, subset=indices,
                    random=list(~1|substudy_num,~1|study_num),method="REML",test="t"),
             silent=TRUE) #with random effects
  
  if (is.element("try-error", class(res))) {
    NA
  } else {
    c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
}

# #for loop to sample by row 
# 
# #fake dataset with test
# t_df <- data.frame("cat"=c(rep("A",6),rep("B",4),rep("C",3)),"yi" = round(rnorm(13,4,2),2), "vi"= round(rnorm(13,.1,0.001),2))
# 
# #empty list to dump results into
# lst <- list()
# 
# sim_mat <- replicate(10, as.character(sample(unique(t_df$cat),replace=TRUE)))
# bootvec<-c()
# 
# for(s in 1:ncol(sim_mat)){
#   vec <- sim_mat[,s]
#   lst<- list()
#   
#   
#   for(i in 1:length(vec)){
#     temp <- filter(t_df,cat==vec[i]) #pull out first value 
#     temp$randstudy<-i
#     lst[[length(lst)+1]] <- temp
#     
#   }
#   
#   temp2 <- tibble(lst)%>%
#     unnest(cols = c(lst))
#   
#   bootvec[s] <-coef((rma.mv(yi=yi, V=vi, data=temp2,
#                             random=list(~1|randstudy),method="REML",test="t")))
#   
# }
# 
# mean(bootvec)

