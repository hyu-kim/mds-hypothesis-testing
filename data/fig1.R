library(ggplot2)
## Panel A
v_lambda = c(0, 1, 2, 4, 10)/10
v_linestyle = c('dotted', 'dashed', 'solid', 'blank', 'blank')
v_linesize = c(0.5, 0.4, 0.2, 0.1, 0.1)
# v_color = c('#f8766d','#7cae00','#00bfc4','black','black')
g <- ggplot()

for(i in 1:length(v_lambda)){
  lambda <- v_lambda[i]
  df_log <- read.csv(paste('result/HyperparameterStudy/sim_1',
                           # ifelse(lambda>0, '50', '200'),
                           '/sim_1-fmds-', sprintf('%.2f',lambda), '-log.csv', sep=''))
  if(lambda>=0.4){
    g <- g + geom_point(data=df_log, aes(x=epoch, y=p_z), size=1.25, shape=c(1,2)[i-3])
  }
  # else{
  g <- g + geom_line(data=df_log, aes(x=epoch, y=p_z), size=v_linesize[i], linetype=v_linestyle[i])
  # }
}

g + 
  scale_y_continuous(breaks = (0:5)/5) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('result/fig1_p.pdf', width=2.15, height=3, units='in')


## Panel B
# import
v_lambda_1 = (0:3)/20 # lambdas for 200-iterated data
v_lambda_2 = c(0,1,2,5,10)/10 # lambdas for 50-iterated data
ratio_v = c(100, 1, 5)
lambda = 0.4
df_log <- read.csv(paste('result/HyperparameterStudy/sim_1',
                         # ifelse(lambda>=0.10, '50', '200'),
                         '/sim_1-fmds-', sprintf('%.2f',lambda), '-log.csv', sep=''))

ggplot(df_log) + 
  # geom_point(aes(x=epoch, y=obj_mds), size=0.75, shape=1) + # iter_50 only
  # geom_line(aes(x=epoch, y=obj_confr*0.02), inherit.aes = FALSE, size=0.25, linetype='dashed') + # iter_50 only
  geom_line(aes(x=epoch, y=obj_mds), size=0.25, linetype='solid') +
  geom_point(aes(x=epoch, y=obj_confr*5), inherit.aes = FALSE, size=0.5, shape=16) +
  scale_y_continuous(sec.axis = sec_axis(~./5, name="Confirmatory")) +
  # ylim(c(0, 6000))+
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave(paste('result/fig1_', sprintf('%.2f', lambda), '_obj.pdf', sep=''), 
       width=2.6, height=1.2, units='in')