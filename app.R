source("modified_grind.R")

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N
    return(list(c(dR, dN)))  
  }) 
}  
p <- c(r=1,K=1,a=1,c=1,delta=0.5) 
s <- c(R=1,N=0.1)

library(shiny)
library(mathjaxr)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-color: #f3f3f3;
        font-family: 'Arial', sans-serif;
      }
      .title-panel {
        text-align: center;
        font-family: 'Arial', sans-serif;
        color: #4CAF50;
        margin-bottom: 20px;
      }
      .sidebar {
        background-color: #fff;
        border-radius: 8px;
        padding: 20px;
        box-shadow: 0px 2px 5px rgba(0,0,0,0.1);
      }
      .slider-label {
        font-weight: bold;
        color: #333;
      }
      .main-panel {
        background-color: #fff;
        border-radius: 8px;
        padding: 20px;
        box-shadow: 0px 2px 5px rgba(0,0,0,0.1);
        margin-top: 20px;
      }
      .mathjax {
        font-size: 18px;
        color: #333;
        margin-bottom: 20px;
      }
      .description {
        margin-top: 20px;
        font-size: 16px;
        color: #555;
      }
      .radio-buttons {
        margin-top: 20px;
        margin-bottom: 20px;
      }
      .numeric-input {
        margin-top: 10px;
      }
      .signature {
        text-align: center;
        font-size: 14px;
        color: #888;
        margin-top: 30px;
      }
    "))
  ),
  
  div(
    class = "title-panel",
    titlePanel("Lotka Volterra Model Simulation")
  ),
  sidebarLayout(
    sidebarPanel(
      class = "sidebar",
      tags$div(
        class = "slider-label",
        sliderInput(inputId=names(s[1]), label=names(s[1]), min=0, max=10, step=0.01, value=as.numeric(s[1]))
      ),
      tags$div(
        class = "slider-label",
        sliderInput(inputId=names(s[2]), label=names(s[2]), min=0, max=10, step=0.01, value=as.numeric(s[2]))
      ),
      tags$div(
        class = "slider-label",
        sliderInput(inputId=names(p[1]), label=names(p[1]), min=0, max=p[1]*10, step=p[1]*0.01, value=p[1])
      ),
      tags$div(
        class = "slider-label",
        sliderInput(inputId=names(p[2]), label=names(p[2]), min=0, max=p[2]*10, step=p[2]*0.01, value=p[2])
      ),
      tags$div(
        class = "slider-label",
        sliderInput(inputId=names(p[3]), label=names(p[3]), min=0, max=p[3]*10, step=p[3]*0.01, value=p[3])
      ),
      tags$div(
        class = "slider-label",
        sliderInput(inputId=names(p[4]), label=names(p[4]), min=0, max=p[4]*10, step=p[4]*0.01, value=p[4])
      ),
      tags$div(
        class = "slider-label",
        sliderInput(inputId=names(p[5]), label=names(p[5]), min=0, max=p[5]*10, step=p[5]*0.01, value=p[5])
      )
    ),
    
    mainPanel(
      class = "main-panel",
      h4("Lotka-Volterra Equations", class = "mathjax"),
      withMathJax(
        "$$ \\frac{dR}{dt} = rR \\left(1 - \\frac{R}{K}\\right) - aRN $$",
        "$$ \\frac{dN}{dt} = caRN - \\delta N $$"
      ),
      plotOutput(outputId="grind"),
      div(class = "radio-buttons",
          radioButtons("radio", "Select Plot Type:", 
                       choices = list("Time plot" = 0, "Nullclines" = 1, "Trajectory" = 2, 
                                      "Phase Portrait" = 3, "Steady State" = 4),
                       selected = 0, inline = TRUE)
      ),
      fluidRow(
        tags$div(class = "numeric-input",
                 column(3, numericInput(inputId="tmax", label="Tmax", value=100))
        ),
        tags$div(class = "numeric-input",
                 column(3, numericInput(inputId="tstep", label="Tstep", value=0.1, step=0.1))
        )
      ),
      textOutput("log"),
      div(class = "description",
          h5("Description:"),
          p("The Lotka-Volterra model describes the dynamics of biological systems in which two species interact, 
             one as a predator and the other as prey. The model consists of a pair of first-order, non-linear, 
             differential equations."),
          p("Variables:"),
          tags$ul(
            tags$li("R: Prey population"),
            tags$li("N: Predator population"),
            tags$li("r: Growth rate of prey"),
            tags$li("K: Carrying capacity of the environment for the prey"),
            tags$li("a: Predation rate coefficient"),
            tags$li("c: Conversion factor of consumed prey to predator growth"),
            tags$li("delta: Death rate of predators")
          )
      ),
      div(class = "signature",
          p("Hope to see you soon at campus!!-- Yiğit Efe Enhoş")
      )
    )
  )
)

server <- function(input, output) {
  output$grind <- renderPlot({
    radiob <- input$radio
    for (i in names(s)) s[i] <- input[[i]]
    for (i in names(p)) p[i] <- input[[i]]
    tmax <- input$tmax
    tstep <- input$tstep
    xmax <- 1.05 * max(p["K"], p["delta"] / (p["c"] * p["a"]), s["R"])
    ymax <- 1.05 * max(p["r"] / p["a"], s["N"])
    output$log <- renderText("")
    
    if (radiob == 0) {
      f <- run(tmax=tmax, tstep=tstep, odes=model, state=s, parms=p)
      output$log <- renderText({paste("Ended in R =", round(f[1], 5), "N =", round(f[2], 5), sep=" ")})
    } else if (radiob <= 2) {
      plane(xmax=xmax, ymax=ymax, odes=model, state=s, parms=p)
      if (radiob == 2) {
        f <- run(tmax=tmax, tstep=tstep, odes=model, state=s, parms=p, traject=TRUE)
        output$log <- renderText({paste("Ended in R =", round(f[1], 5), "N =", round(f[2], 5), sep=" ")})
      }
    } else if (radiob == 3) {
      plane(xmax=xmax, ymax=ymax, tmax=tmax, tstep=tstep, odes=model, state=s, parms=p, portrait=TRUE)
    } else {
      plane(xmax=xmax, ymax=ymax, odes=model, state=s, parms=p)
      f <- newton(state=s, parms=p, odes=model, plot=TRUE, positive=TRUE)
      if (is.null(f)) output$log <- renderText({"No convergence: start closer to a steady state"})
      else output$log <- renderText({paste("Converged into R =", round(f[1], 5), "N =", round(f[2], 5), sep=" ")})
    }
  })
}

shinyApp(ui = ui, server = server)
