trans_irf = function(dat){
  dat = dat %>% sapply(function(x) x) %>% t %>% as.data.frame() 
  dat$V1 = cumsum(dat$V1)
  dat$V4 = cumsum(dat$V4)
  dat$V7 = cumsum(dat$V7)
  dat
}

require(ggplot2)
plot_irf = function(dat, variable = c('oil_prod', 'real_act', 'oil_price')){
  
  h = 0:(nrow(dat)-1)
  dat = data.frame(h, dat) %>% reshape2::melt('h')
  dat$variable = factor(dat$variable, levels = c('V1', 'V4', 'V7', 'V2', 'V5', 'V8', 'V3', 'V6', 'V9'))
  
  my_labeller <- as_labeller(c(V1 = paste("os", "%->%", variable[1]), 
                               V2 = paste("os", "%->%", variable[2]), 
                               V3 = paste("os", "%->%", variable[3]), 
                               V4 = paste("ad", "%->%", variable[1]), 
                               V5 = paste("ad", "%->%", variable[2]), 
                               V6 = paste("ad", "%->%", variable[3]), 
                               V7 = paste("osd", "%->%", variable[1]), 
                               V8 = paste("osd", "%->%", variable[2]), 
                               V9 = paste("osd", "%->%", variable[3])), 
                             default = label_parsed)
  
  dat %>% ggplot(aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
    facet_wrap(~variable, scales = "free_y", ncol = 3, labeller = my_labeller) + xlab(" ") + ylab(" ") +
    theme_bw() 
}
