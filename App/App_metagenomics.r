#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

list.of.packages <- c("shiny","shinycssloaders","phyloseq","vegan","rstatix","ggplot2", "pscl")
new.packages <- setdiff(list.of.packages,installed.packages()[,"Package"])
if(!require(BiocManager)){
  install.packages("BiocManager")
  library(BiocManager)
}

if(length(new.packages)>0){ install(new.packages,update=TRUE,ask=FALSE)}

library(shiny)
library(shinycssloaders)
library(phyloseq)
library(vegan)
library(rstatix)
library(ggplot2)
library(pscl)

ui <- navbarPage("Metagenomic analysis",
                 tabPanel("Load and apply general or specific settings",
                          tabsetPanel(
                            tabPanel("Load your data",
                                     fluidRow(
                                       column(12,
                                              fileInput("file",label = "Select or drop your phyloseq object",accept=".rds")),
                                       column(3,
                                              uiOutput("select_variable")),
                                       column(3,
                                              uiOutput("select_group")),
                                       column(3,
                                              uiOutput("select_rank"),)
                                     )
                            ),
                            # Panel for alpha-diversity
                            tabPanel("Alpha diversity settings",
                                     selectInput("alpha","Alpha diversity index",c("Chao1","Shannon","Simpson"))
                                         ),
                            # Panel for beta-diversity
                            tabPanel("Beta diversity settings",
                                     selectInput("transpose",label = "Are the taxa as rows in your otu table?",c("Yes","No"),selected = "No"),
                                     selectInput("beta",label = "Beta diversity method",c("bray", "jaccard", "aitchison","weighted unifrac","unweighted unifrac")),
                                     sliderInput("percentage_beta",
                                                 label = "Percentage of prevalence",
                                                 min = 0,
                                                 max = 100,
                                                 value = 0)
                                     ),
                            # Panel for ZINB
                            tabPanel("ZINB settings",
                                     sliderInput("percentage_zinb",
                                                 label = "percentage of prevalence in all samples",
                                                 min = 0,
                                                 max = 100,
                                                 value = 0),
                                     sliderInput("percentage_zinb_group",
                                                 label = "percentage of prevalence in group samples",
                                                 min = 0,
                                                 max = 100,
                                                 value = 0),
                                     textInput("qvalue",
                                                 label = "q-value cut off for plotting",
                                                  "0.05"),
                                     textInput("foldchange",
                                               label = "Foldchange for plotting","1.5"))
                          )
                 ),
                 tabPanel("Alpha diversity",
                          tabsetPanel(
                            tabPanel("Plot",
                                     withSpinner(plotOutput(outputId = "alpha_plot"),color = getOption("spinner.color", default = "dodgerblue3"),
                                                 type=getOption("spinner.type",default=2),color.background = "dodgerblue4")),
                            tabPanel("Table",
                                     h2("Results from kruskal walis"),
                                     withSpinner(dataTableOutput(outputId = "table_k_alpha"),
                                                 color = getOption("spinner.color", default = "dodgerblue3"),
                                                 type=getOption("spinner.type",default=2),color.background = "dodgerblue4"),
                                     h2("Results from glm gamma"),
                                     withSpinner(dataTableOutput(outputId = "table_glm_alpha"),color = getOption("spinner.color", default = "dodgerblue3"),
                                                 type=getOption("spinner.type",default=2),color.background = "dodgerblue4")
                          )
                 )),
                 tabPanel("Beta diversity",
                          tabsetPanel(
                            tabPanel("Plot",
                                     withSpinner(plotOutput(outputId = "beta_plot"),
                                                 color = getOption("spinner.color", default = "dodgerblue3"),
                                                 type=getOption("spinner.type",default=2),color.background = "dodgerblue4")),
                            tabPanel("Table",
                                     withSpinner(dataTableOutput("beta_table"),
                                                 color = getOption("spinner.color", default = "dodgerblue3"),
                                                 type=getOption("spinner.type",default=2),color.background = "dodgerblue4"))
                          )
                          ),
                 tabPanel("ZINB",
                          tabsetPanel(
                            tabPanel("Plot",
                                     withSpinner(plotOutput(outputId = "ZINB_plot"),
                                                 color = getOption("spinner.color", default = "dodgerblue3"),
                                     type=getOption("spinner.type",default=2),color.background = "dodgerblue4")),
                            tabPanel("Table",
                                     withSpinner(dataTableOutput(outputId = "ZINB_table"),color = getOption("spinner.color", default = "dodgerblue3"),
                                                 type=getOption("spinner.type",default=2),color.background = "dodgerblue4")))
                          )
                 )

  
server <- function(input, output,session) {
observeEvent(input$file,{
  if (!is.null(input$file)) {
    archivo<-readRDS(input$file$datapath)
    clinica<-data.frame(sample_data(archivo))
    output$select_variable<-renderUI({
      selectInput("variable","Select your grouping_variable",colnames(clinica))
    })
  }
})
  
  observeEvent(input$variable,{ 
    if (!is.null(input$variable)) {
      archivo<-readRDS(input$file$datapath)
      clinica<-data.frame(sample_data(archivo))
      output$select_group<-renderUI({
        selectInput("groups","Select your groups",unique(clinica[,input$variable]),multiple=TRUE)
      })
    }
  })
  
  observeEvent(input$file,{
    if (!is.null(input$file)) {
      archivo<-readRDS(input$file$datapath)
      taxa<-data.frame(tax_table(archivo))
      output$select_rank<-renderUI({
        selectInput("rank","Select taxonomic rank for beta and relative abundance analysis",colnames(taxa))
      })
    }
  })

  # Alpha diversity
  output$alpha_plot<-renderPlot({
    req(input$groups,input$variable)
    ps<-readRDS(input$file$datapath)
    alpha_div_ps<-estimate_richness(ps,measures=c("Chao1","Shannon","Simpson"))
    alpha_div_ps$grouping_variable<-as.factor(data.frame(sample_data(ps))[,input$variable])
    alpha_div_ps_mod<-alpha_div_ps[alpha_div_ps$grouping_variable==input$groups,]
    
    
    ggplot(alpha_div_ps_mod,aes(x=grouping_variable,y=get(input$alpha)))  +
      geom_violin(aes(colour=grouping_variable)) + geom_jitter(width = 0.2,aes(colour=grouping_variable)) +
      labs(y=input$alpha) +
      scale_fill_discrete(name=input$variable) + theme_minimal()
  })
  output$table_k_alpha<-renderDataTable({
    req(input$groups,input$variable)
    ps<-readRDS(input$file$datapath)
    alpha_div_ps<-estimate_richness(ps,measures=c("Chao1","Shannon","Simpson"))
    alpha_div_ps$grouping_variable<-as.factor(data.frame(sample_data(ps))[,input$variable])
    alpha_div_ps_mod<-alpha_div_ps[alpha_div_ps$grouping_variable%in%input$groups,]
    
    dunn_test(alpha_div_ps[alpha_div_ps$grouping_variable%in%input$groups,],as.formula(paste0(input$alpha,"~","grouping_variable")),p.adjust.method = "fdr")
  })
  output$table_glm_alpha<-renderDataTable({
    req(input$groups,input$variable)
    ps<-readRDS(input$file$datapath)
    alpha_div_ps<-estimate_richness(ps,measures=c("Chao1","Shannon","Simpson"))
    alpha_div_ps$grouping_variable<-as.factor(data.frame(sample_data(ps))[,input$variable])
    alpha_div_ps_mod<-alpha_div_ps[alpha_div_ps$grouping_variable==input$groups,]
    glm_res<-summary(glm(as.formula(paste0(input$alpha,"~","grouping_variable")),data=alpha_div_ps[alpha_div_ps$grouping_variable%in%input$groups,],family = Gamma(link="log")))$coefficients
    glm_res<-cbind(rownames(glm_res),glm_res)
    as.data.frame(glm_res)
    })
  
  # Beta diversity
datos_beta<-reactive({
  req(length(input$groups)>1,input$variable)
    ps<-readRDS(input$file$datapath)
    if (input$transpose=="Yes") {
      otu_table(ps)<-t(otu_table(ps))
    }
    min_samples <- round((nrow(sample_data(ps))*as.numeric(input$percentage_beta)/100)) 
    ps<-tax_glom(ps,taxrank = input$rank)
    data_phylo_filt <- prune_taxa(taxa_sums(ps) >= min_samples, ps)
    sample_data(data_phylo_filt)$grouping_variable<-as.factor(data.frame(sample_data(data_phylo_filt))[,input$variable])
    taxa<-tax_table(data_phylo_filt)
    samples_df<-sample_data(data_phylo_filt)
    set.seed(1782) # set seed for analysis reproducibility
    OTU_filt_rar = rarefy_even_depth(otu_table(data_phylo_filt), rngseed = TRUE, replace = FALSE) # rarefy the raw data using Phyloseq package
    data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar)) # create a separated file
    data_phylo_filt_rar <- phyloseq(OTU_filt_rar, taxa, sample_data(samples_df)) # create a phyloseq object
    
    ps_temp<-phyloseq(otu_table(data_phylo_filt_rar)[sample_data(data_phylo_filt_rar)$grouping_variable %in% input$groups,],
                      taxa,
                      sample_data(data_phylo_filt_rar)[sample_data(data_phylo_filt_rar)$grouping_variable %in% input$groups,])
    ps_clr<-microbiome::transform(data_phylo_filt,"clr")
    list(ps_temp=ps_temp,
         ps_clr=ps_clr,
         data_phylo_filt=data_phylo_filt)
})
  
  output$beta_plot<-renderPlot({
    req(datos_beta())
    
    ps_temp<-datos_beta()$ps_temp
    ps_clr<-datos_beta()$ps_clr
    data_phylo_filt<-datos_beta()$data_phylo_filt
    if (input$beta=="aitchison") {
      
      
      
      ps_clr_temp<-phyloseq(otu_table(ps_clr)[sample_data(ps_clr)$grouping_variable %in% input$groups,],
                            taxa,
                            sample_data(ps_clr)[sample_data(ps_clr)$grouping_variable %in% input$groups,])
      
      dist_matrix <- vegdist(otu_table(ps_clr_temp), method = "euclidean")
      
      
      # PCoA)
      ordination_result <- pcoa(dist_matrix)
      
      # plot_ordination
      plot_ordination(ps_clr_temp, ordination_result, color = "grouping_variable") + geom_point(size=3) + stat_ellipse() + theme_minimal()
      
    }
    else if (input$beta=="weighted unifrac") {
      tree<-phy_tree(data_phylo_filt)
      phy_tree(ps_temp)<-tree
      
      ordu.wt.uni <- ordinate(ps_temp, "PCoA", "unifrac", weighted=T)
      
      
      wt.unifrac <- plot_ordination(ps_temp, 
                                    ordu.wt.uni, color="grouping_variable") 
      wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
      wt.unifrac <- wt.unifrac + theme_classic() + scale_color_brewer("grouping_variable", palette = "Set2")
      print(wt.unifrac+ stat_ellipse()) + theme_minimal()
      
      
      unifrac.dist <- UniFrac(ps_temp, 
                              weighted = TRUE, 
                              normalized = TRUE,  
                              parallel = FALSE, 
                              fast = TRUE)
      
      permanova <- adonis2(unifrac.dist ~ sample_data(ps_temp)$grouping_variable,na.action = na.omit)
      
      permanova
      
      a<-betadisper(unifrac.dist, sample_data(ps_temp)$grouping_variable, type = c("median","centroid"), bias.adjust = FALSE,
                    sqrt.dist = FALSE, add = FALSE)
      print(as.data.frame(a$grouping_variable.distances))
    }
    
    else if (input$beta=="unweighted unifrac") {
      tree<-phy_tree(data_phylo_filt)
      phy_tree(ps_temp)<-tree
      
      ordu.unwt.uni <- ordinate(ps_temp, "PCoA", "unifrac", weighted=F)
      unwt.unifrac <- plot_ordination(ps_temp, 
                                      ordu.unwt.uni, color="grouping_variable") 
      unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
      unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("grouping_variable", palette = "Set2")
      print(unwt.unifrac+ stat_ellipse()) + theme_minimal()
      
      unifrac.dist <- UniFrac(ps_temp, 
                              weighted = FALSE, 
                              normalized = TRUE,  
                              parallel = FALSE, 
                              fast = TRUE)
      
      permanova <- adonis2(unifrac.dist ~ sample_data(ps_temp)$grouping_variable,na.action = na.omit)
      
      permanova
      
      a<-betadisper(unifrac.dist, sample_data(ps_temp)$grouping_variable, type = c("median","centroid"), bias.adjust = FALSE,
                    sqrt.dist = FALSE, add = FALSE)
      print(as.data.frame(a$grouping_variable.distances))
      
    }
    else{
      pcoa_bc_temp=ordinate(ps_temp,"PCoA",input$beta)
      plot_ordination(ps_temp, pcoa_bc_temp, color = "grouping_variable") +
        geom_point(size = 3) + stat_ellipse() + theme_minimal()
      
    }
  })
  
  output$beta_table<-renderDataTable({
    req(input$groups,input$variable)
    ps_temp<-datos_beta()$ps_temp
    ps_clr<-datos_beta()$ps_clr
    data_phylo_filt<-datos_beta()$data_phylo_filt
    
    if (input$beta=="aitchison") {
      a<-as.data.frame(
        adonis2(as.formula(otu_table(ps_clr)[sample_data(ps_clr)$grouping_variable %in% input$groups,]~(sample_data(ps_clr))$grouping_variable[sample_data(ps_clr)$grouping_variable%in%input$groups]),permutations=9999, method="euclidean",na.action = na.omit)
        
      )
      rownames(a)[1]<-as.character(input$variable)
      a
    }
    else if (input$beta=="weighted unifrac") {
      tree<-phy_tree(data_phylo_filt)
      phy_tree(ps_temp)<-tree
      
      unifrac.dist <- UniFrac(ps_temp,
                              weighted = TRUE,
                              normalized = TRUE,
                              parallel = FALSE,
                              fast = TRUE)
      
      permanova <- adonis2(unifrac.dist ~ sample_data(ps_temp)$grouping_variable,na.action = na.omit)
      
      as.data.frame(permanova)
      
    }
    else if (input$beta=="unweighted unifrac") {
      tree<-phy_tree(data_phylo_filt)
      phy_tree(ps_temp)<-tree
      
      unifrac.dist <- UniFrac(ps_temp,
                              weighted = FALSE,
                              normalized = TRUE,
                              parallel = FALSE,
                              fast = TRUE)
      
      permanova <- adonis2(unifrac.dist ~ sample_data(ps_temp)$grouping_variable,na.action = na.omit)
      
      as.data.frame(permanova)
    }
    else{
      a<-       data.frame(
        adonis2(as.formula(otu_table(ps_temp)[sample_data(ps_temp)$grouping_variable %in% input$groups,]~(sample_data(ps_temp))$grouping_variable[sample_data(ps_temp)$grouping_variable%in%input$groups]),permutations=9999, method=input$beta,na.action = na.omit)
      )
      rownames(a)[1]<-as.character(input$variable)
      a
    }
  })
# ZINB
  datos_zinb<-reactive({
    req(input$groups,input$variable)
    ps<-readRDS(input$file$datapath)
    if(input$transpose=="Yes"){
      otu_table(ps)<-t(otu_table(ps))
    }
    min_samples <- round((nrow(sample_data(ps))*as.numeric(input$percentage_zinb)/100))
    sample_data(ps)$total.counts<-rowSums(otu_table(ps))
    ps<-tax_glom(ps,taxrank = input$rank)
    data_phylo_filt <- prune_taxa(taxa_sums(ps) >= min_samples, ps)
    sample_data(data_phylo_filt)$grouping_variable<-as.factor(data.frame(sample_data(data_phylo_filt))[,input$variable])
    taxa<-tax_table(data_phylo_filt)
    samples_df<-sample_data(data_phylo_filt)
    set.seed(1782)
    data_otu_filt_rar = data.frame(otu_table(data_phylo_filt)) # create a separated file
    data_phylo_filt_rar <- phyloseq(otu_table(data_phylo_filt), taxa, sample_data(samples_df)) # create a phyloseq object
    
    ps_mod2<-phyloseq(otu_table(data_phylo_filt_rar)[sample_data(data_phylo_filt_rar)$grouping_variable %in% input$groups,],
                      taxa,
                      sample_data(data_phylo_filt_rar)[sample_data(data_phylo_filt_rar)$grouping_variable %in% input$groups,])
    min_samples_group <- round((nrow(sample_data(ps_mod2))*as.numeric(input$percentage_zinb_group)/100))  # Ajusta según tu        
    ps_mod2 <- prune_taxa(taxa_sums(ps_mod2) >= min_samples_group, ps_mod2)
    taxa<-tax_table(ps_mod2)
    N=as.numeric(sample_data(ps_mod2)$total.counts)
    pheno=round(as.data.frame(otu_table(ps_mod2))/100*N)
    yy = as.matrix(pheno)  
    yy = ifelse(is.na(yy), 0, yy)
    zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
    zero.p = sort(zero.p, decreasing = T)
    zero.p = data.frame(zero.p)
    zero.p$id = rownames(zero.p)
    zero.p = data.frame(zero.p[zero.p$zero.p>0.2 & zero.p$zero.p<0.8, ])
    yy = yy[, rownames(zero.p)]
    taxa2<-taxa[rownames(zero.p),]
    
    ps_zinb <- phyloseq(otu_table(yy,taxa_are_rows = FALSE), taxa2, sample_data(ps_mod2)) # create a phyloseq object
    adjust<-offset(log(N))
    datos<-cbind(sample_data(ps_zinb),adjust,yy)
    colnames(datos)<-c(colnames(sample_data(ps_zinb)),"adjust",taxa2[,input$rank])
    colnames(datos)<-gsub("-","_",colnames(datos))
    datos<-as.data.frame(datos)
    taxa_glm<-gsub("-","_",taxa2[,input$rank])
    tabla<-NULL
    for (i in taxa_glm) {
      m1<- zeroinfl(datos[,i] ~ sample_data(datos)$grouping_variable ,
                    dist = "negbin")
      res<-summary(m1)
      p<-res$coefficients$zero[2,4]
      estimate<-round(exp(res$coefficients$zero[2,1]),2)
      IC2.5<-round(exp(confint(m1)[4,1]),2)
      IC97.5<-round(exp(confint(m1)[4,2]),2)
      tabla<-rbind(tabla,c(estimate,IC2.5,IC97.5,p))
    }
    rownames(tabla)<-taxa2[,input$rank]
    
    tabla<-as.data.frame(tabla)
    colnames(tabla)<-c("estimate","IC2.5","IC97.5","p")
    tabla$p.adj<-p.adjust(as.numeric(tabla$p),method = "fdr")
    tabla$p.adj<-round(tabla$p.adj,3)
    tabla$p<-round(tabla$p,3)
    tabla<-cbind(rownames(tabla),tabla,as.data.frame(tax_table(ps_zinb))$Phylum)
    colnames(tabla)<-c(input$rank,"Estimate","IC 2.5","IC97.5","p","p.adj (fdr)","Phylum")
    tabla<-as.data.frame(tabla)
    tabla2<-tabla[tabla$"p.adj (fdr)"<=as.numeric(input$qvalue),]
    
      p<-ggplot(tabla2, aes(y=get(input$rank), x=log(Estimate), color=Phylum)) +
        geom_vline(xintercept = c(-log2(as.numeric(input$foldchange)),0,log2(as.numeric(input$foldchange))), color = c("red","gray","red"), size = 0.5) +
        geom_point(size=6) +
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +theme_minimal()

    list(table_zinb=tabla,
         graph=p)
  })
  output$ZINB_plot<-renderPlot({
    p<-datos_zinb()$graph
    p
  })
  output$ZINB_table<-renderDataTable({
    table_zinb<-data.frame(datos_zinb()$table_zinb)
    table_zinb
  })
  
  
  }

# Run the application 
shinyApp(ui = ui, server = server)

