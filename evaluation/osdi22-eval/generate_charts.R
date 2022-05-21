#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

args <- c('all-data-alt-names.csv')
args <- commandArgs(trailingOnly=TRUE)
bench <- data.frame(read.csv(args[1], header=FALSE, 
                             col.names = c('suite','script','shell','speedup')))

text <- element_text(family='Times', size=17)

box_order <- c('Classics', 'Unix50', 'COVID-mts', 'NLP')
boxes <- 
  ggplot(bench[bench$shell %in% c('pash_aot', 'pash_jit') & !(bench$suite %in% c('for-loops', 'AvgTemp', 'WebIndex')),], 
         aes(x=factor(suite, level=box_order), y=as.double(speedup), fill=factor(rev(shell)))) +
  geom_boxplot() + #geom(aes(color=factor(shell))) +
  geom_hline(yintercept = 1) +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(family='Times'),
        legend.position = 'none') +
  labs(y = "Speedup vs. bash (log on left)", x = "", color = "Shell")
bar_order <- c('AvgTemp', 'WebIndex', 'MediaConv1', 'MediaConv2', 'ProgInf', 'LogAnalysis1', 'LogAnalysis2', 'Genomics', 'AurPkg', 'FileEnc1', 'FileEnc2')
bars <- 
  ggplot(bench[bench$shell %in% c('pash_aot', 'pash_jit') & bench$suite %in% c('for-loops', 'AvgTemp', 'WebIndex'),],
         aes(x=factor(script, level=bar_order), y=as.double(speedup), fill=factor(rev(shell)))) +
  geom_col(position=position_dodge2(reverse=TRUE)) +
  geom_hline(yintercept = 1) +
  labs(y = "", x = "", fill = "") +
  scale_fill_discrete(labels = c('PaSh JIT', 'PaSh AOT')) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(family='Times'),
        legend.position = c(0.9, 0.81),
        legend.title = element_blank()) + 
      guides(fill = guide_legend(reverse = FALSE))
p <- grid.arrange(boxes, bars, nrow=1, ncol=2)
ggsave("figure5.pdf", p, width=11, height=3)

boxes <-
  ggplot(bench[bench$shell %in% c('pash_jit', 'pash_jit -prof -par_pipe', 'pash_jit -prof') & bench$suite %in% c('NLP'),],
         aes(x=factor(suite), y=as.double(speedup), fill=factor(shell))) +
  geom_boxplot(aes(x=factor(suite, level=box_order)), position=position_dodge2(reverse=TRUE), key_glyph = draw_key_rect) +
  #geom_point(aes(x=factor(suite, level=box_order), color=factor(shell)), position = position_dodge2(width=0.75, reverse=TRUE), show.legend = FALSE) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(labels = c('NLP')) +
  labs(y = "Speedup vs. bash", x = "", fill="") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = text,
        legend.position='none')
bars <-
  ggplot(bench[bench$shell %in% c('pash_jit', 'pash_jit -prof -par_pipe', 'pash_jit -prof') & bench$suite %in% c('for-loops','AvgTemp'),],
         aes(x=factor(script, level=bar_order), y=as.double(speedup), fill=factor(shell))) +
  geom_col(position=position_dodge2(reverse=TRUE, width=1)) +
  geom_hline(yintercept = 1) +
  labs(y = "", x = "", fill = "") +
  scale_fill_discrete(labels = c('PaSh JIT', 'PaSh JIT no_prof', 'PaSh JIT no_prof no_du')) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = text,
        legend.position = c(0.79, 0.80),
        legend.title = element_blank())
p <- grid.arrange(boxes, bars, nrow=1, ncol=2, widths = c(1,4))
ggsave("figure6.pdf", p, width=8, height=4)



# data munging to get the right names
# bench$cshell <- bench$shell
# bench[!(bench$cshell %in% c('pash_jit', 'pash_jit_no_comm')),]$cshell <- 'ignore'
# bench[bench$shell == 'pash_jit_no_comm',]$cshell <- 'pash_jit no_comm'
# bench[bench$shell == 'pash_jit',]$cshell <- 'pash_jit'

# bench$shell %in% c('pash_jit', 'pash_jit -prof -par_pipe', 'pash_jit -prof')

boxes <- 
  ggplot(bench[bench$shell %in% c('pash_jit', 'pash_jit_no_comm') & bench$suite %in% c('Classics','Unix50','COVID-mts'),],
         aes(y=as.double(speedup), fill=factor(shell))) +
#   ggplot(bench[bench$cshell != 'ignore' & bench$suite %in% c('Classics','Unix50','COVID-mts'),],
#          aes(y=as.double(speedup), fill=factor(cshell))) +
  geom_boxplot(aes(x=factor(suite, level=box_order)), position=position_dodge2(reverse=TRUE)) +
  #geom_point(aes(x=factor(suite, level=box_order), color=factor(cshell)), position = position_dodge2(width=0.75, reverse=TRUE)) +
  geom_hline(yintercept = 1) +
  labs(y = "Speedup vs. bash", x = "", fill="", color="") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = text,
        legend.position = 'none')
bars <- 
  ggplot(bench[bench$shell %in% c('pash_jit', 'pash_jit_no_comm') & bench$suite %in% c('AvgTemp','WebIndex'),],
         aes(x=factor(script, level=bar_order), y=as.double(speedup), fill=factor(shell))) +
#   ggplot(bench[bench$cshell != 'ignore' & bench$suite %in% c('AvgTemp','WebIndex'),],
#          aes(x=factor(script, level=bar_order), y=as.double(speedup), fill=factor(cshell))) +
  geom_col(position=position_dodge2(reverse=TRUE)) +
  geom_hline(yintercept = 1) +
  scale_fill_discrete(labels = c('PaSh JIT', 'PaSh JIT no_comm')) +
  labs(y = "", x = "", fill="") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = text,
        legend.position = c(0.27, 0.82),
        legend.title = element_blank())
p <- grid.arrange(boxes, bars, nrow=1, ncol=2)
ggsave("figure7.pdf", p, width=8, height=4)

