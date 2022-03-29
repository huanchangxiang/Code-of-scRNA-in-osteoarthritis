
PropPlot <- function(object, groupBy){
  # (1)获取绘图数据
  plot_data = object@meta.data %>% 
    dplyr::select(orig.ident, {{groupBy}}) %>% 
    dplyr::rename(group = as.name(groupBy))
  
  # (2)绘图
  figure = ggstatsplot::ggbarstats(data = plot_data, 
                      x = group, y = orig.ident,
                      package = 'ggsci',
                      palette = 'category20c_d3',
                      results.subtitle = FALSE,
                      bf.message = FALSE,
                      proportion.test = FALSE,
                      label.args = list(size = 5, 
                                        fill = 'white', 
                                        alpha = 0.85,
                                        family = 'Arial',
                                        fontface = 'bold'),
                      perc.k = 2,
                      title = '',
                      xlab = '',
                      legend.title = "Cluster",
                      ggtheme = ggpubr::theme_pubclean()) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = 'black', lineend = 'round'),
          legend.position = 'right',
          axis.text.x = element_text(size = 15, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 15, color = 'black', family = 'Arial'),
          legend.text = element_text(family = 'Arial', size = 15, color = 'black'),
          legend.title = element_text(family = 'Arial', size = 17, color = 'black')) 
  
  # (3)去除柱子下面的样本量标识：
  gginnards::delete_layers(x = figure, match_type = 'GeomText')
}


PropPlot(object = experiment.aggregate, groupBy = 'RNA_snn_res.0.5')+scale_fill_manual(values = c("#CC99FFFF","#00FFFFFF","#95CC5EFF","#CC0C00FF"))
