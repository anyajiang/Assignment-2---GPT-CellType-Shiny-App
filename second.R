library(shiny)
library(Seurat)
library(openai)

Sys.setenv(OPENAI_API_KEY = 'XXXXX')

# Example data
data("pbmc_small")

# UI for application 
ui <- fluidPage(
  titlePanel("GPTCelltype Annotation"), 
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select Seurat Differential Gene Table", 
                  choices = list("pbmc_small" = "1", "input list" = "2", "others" = "3"), 
                  selected = "pbmc_small"),
      selectizeInput("model", label = "Select Model", 
                     choices = c("gpt-3.5-turbo", "gpt-4"), 
                     selected = "gpt-4", multiple = FALSE),
      conditionalPanel(
        condition = "input.dataset == '2'",
        textAreaInput("gene_list", "Enter Gene List (New line separated for each group. Space separated for each gene)", 
                      "", rows = 5), 
        actionButton("example", "Example Query")
      ),
      actionButton("annotate", "Annotate Cell Types")
    ), 
    mainPanel(
      titlePanel("Cell Type Annotation: "), 
      conditionalPanel(
        condition = "input.dataset == '1'", 
        plotOutput("dimPlot")
      ),
      conditionalPanel(
        condition = "input.dataset == '2'", 
        textOutput("message")
      ), 
      conditionalPanel(
        condition = "input.dataset != '1' & input.dataset != '2'", 
        textOutput("other")
      )
    )
  )
)

# Server for application 
server <- function(input, output, session) {
  
  observeEvent(input$annotate, {
    output$message <- NULL
    output$dimPlot <- NULL 
    output$other <- NULL
    
    if (input$dataset == 1) {
      dataset <- pbmc_small
      
      suppressWarnings({
        all.markers <- FindAllMarkers(object = dataset)
      })
      res <- gptcelltype(all.markers, 
                         tissuename = 'human PBMC', 
                         model = input$model)
      dataset@meta.data$celltype <- as.factor(res[as.character(Idents(dataset))])
      
      output$dimPlot <- renderPlot({
        DimPlot(dataset, group.by = 'celltype')
      })
    } else if (input$dataset == 2) {
      gene_groups <- strsplit(input$gene_list, "\n")[[1]]
      gene_list <- lapply(gene_groups, function(x) strsplit(x, ",")[[1]])
      names(gene_list) <- paste0("Group", seq_along(gene_list))
      
      res <- gptcelltype(input = gene_list, tissuename = NULL, model = input$model, topgenenumber = 10) 
      output$message <- renderPrint(
        res
      )
    } else {
      output$other <- renderTet({
        "We only support the first two options!"
      })
    }
  })
  
  observeEvent(input$example, {
    updateTextAreaInput(session, "gene_list", value = "CD4 CD3D\nCD14")
  })
}

# GPT cell type function 
gptcelltype <- function(input, tissuename = NULL, model = 'gpt-4', topgenenumber = 10) {
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  if (OPENAI_API_KEY == "") {
    print("Note: OpenAI API key not found: returning the prompt itself.")
    API.flag <- 0
  } else {
    API.flag <- 1
  }
  
  if (class(input) == 'list') {
    input <- sapply(input, paste, collapse = ',')
  } else {
    input <- input[input$avg_log2FC > 0,, drop = FALSE]
    input <- tapply(input$gene, list(input$cluster), function(i) paste0(i[1:topgenenumber], collapse = ','))
  }
  
  if (!API.flag) {
    message <- paste0('Identify cell types of ', tissuename, ' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types. ',  "\n", paste0(names(input), ':', unlist(input), collapse = "\n"))
    return(message)
  } else {
    print("Note: OpenAI API key found: returning the cell type annotations.")
    cutnum <- ceiling(length(input) / 30)
    if (cutnum > 1) {
      cid <- as.numeric(cut(1:length(input), cutnum))	
    } else {
      cid <- rep(1, length(input))
    }
    
    allres <- sapply(1:cutnum, function(i) {
      id <- which(cid == i)
      flag <- 0
      while (flag == 0) {
        k <- openai::create_chat_completion(
          model = model,
          message = list(list("role" = "user", "content" = paste0('Identify cell types of ', tissuename, ' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types.\n', paste(input[id], collapse = '\n'))))
        )
        res <- strsplit(k$choices[,'message.content'], '\n')[[1]]
        if (length(res) == length(id))
          flag <- 1
      }
      names(res) <- names(input)[id]
      res
    }, simplify = F) 
    print('Note: It is always recommended to check the results returned by GPT-4 in case of\n AI hallucination, before going to down-stream analysis.')
    return(gsub(',$', '', unlist(allres)))
  }
}

# Run the app 
shinyApp(ui, server)
