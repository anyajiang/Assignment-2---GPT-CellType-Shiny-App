library(shiny)
library(openai)
#library(jsonlite)
Sys.setenv(openai_api_key = 'sk-proj-lbwdwtz8NFqxYN9Uq0mzT3BlbkFJ5ZaMHWFZIEruufHN68ry')  #Enter your openAI API key here
chatHist <- t(list())
AI_content <- reactiveVal(list())
# Define UI for application 
ui <- fluidPage(
  tags$head(
    tags$style(type="text/css", "#assistant{ height: 800px; font-family: sans-serif; white-space: pre-wrap;}", 
              ".shiny-notification{
                 position: fixed;
                 top: 33%;
                 left: 33%;
                 right: 33%;
               }")
  ),
  titlePanel("Chat Interface"),
  sidebarLayout(
    sidebarPanel(
      textAreaInput("user_message", "Enter your message:", rows=6),
      sliderInput("temperature", "Temperature (low: focused, high: creative):", min = 0, max = 2, value = 0.7, step = 0.1),
      selectizeInput("model", label="Select Model", 
                     choices=c(c("gpt-3.5-turbo", "gpt-4", "gpt-4-turbo-preview", "gpt-4o")), selected="gpt-4o", multiple=F),
      actionButton("submit_message", "Submit Message")
    ),
    mainPanel(
      # in the UI function
      titlePanel("Response:"),
      span(textOutput("chat_rounds"), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
      downloadButton("download_list", "Download Chat History"),
      actionButton("clear_history", "Clear Chat History"),
      textOutput(outputId = "assistant")    
    )
  )
)

server <- function(input, output, session) {
  
  # Store chat history
  observeEvent(input$clear_history, {
    chatHist(list())
    AI_content(list())
    updateTextInput(session, "user_message", value = "")
  })
  
  observeEvent(input$submit_message, {
    user_message <- input$user_message
    chat_history <- chatHist() 
    # Add user's message to chat history
    chat_history <- append(chat_history, list(list("role" = "user",
                                                   "content" = user_message)))
    
    # Clear the user input field
    updateTextInput(session, "user_message", value = "")
    
    # Call the create_chat_completion function
    withProgress(message = 'Generating response using OpenAI...', value = 0, {
    assistant_response <- create_chat_completion(
      model = input$model,
      messages = chat_history,
      temperature = input$temperature,
      openai_api_key = Sys.getenv("openai_api_key")
    ) 
    })
    # Extract and append assistant's response content to chat history
    assistant_content <- assistant_response$choices
    chat_history <- append(chat_history, list(list("role" = "assistant", 
                                                   "content" = assistant_content$message.content)))
    chatHist(chat_history)
    assistant_response$user_message=user_message
    AI_content(assistant_response)
    # Update the chat display in the mainPanel
    #browser()
  })
  
  observe({
    chat_history <- chatHist() 
    assistant_response=AI_content()
    output$chat_rounds <- renderText({paste("Conversation length: ",  ceiling(length(chat_history)/2), " rounds", sep="")})
    assistant_content <- assistant_response$choices
    output$assistant <-renderText({
        if (length(assistant_response)==0) {
          "\n\nType a message and click submit to start"
        } else {
          paste0("\n", assistant_content$message.content, "\n\nYour Prompt: ",assistant_response$user_message,  "\n\nModel: ",  assistant_response$model, "\n\nToken use: ", assistant_response$usage$total_tokens )
        }
    })
  })
  
  output$download_list <- downloadHandler(
    filename = function() {
      "chat_history.txt"
    },
    content = function(file) {
      chat_history <- chatHist() 
      for (i in 1:length(chat_history)){
        L1=chat_history[[i]]
        cat("##", L1[[1]], "\n", L1[[2]], "\n\n", sep="", file=file, append=T)
      }
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
