library(shiny)
library(ggplot2)
library(splines)

load("dp_ranef_s.RData")

# Knots placement
inner_knots <- 40; degree <- 3
xl    <- 15; xr <- 49; dx <- (xr - xl) / (inner_knots-1)
knots <- seq(xl - degree * dx, xr + degree * dx, by = dx)

# Fixed quantities
B     <- spline.des(knots, 15:49, degree + 1, 0 * 15:49, outer.ok=TRUE)$design


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  prob1 <- reactive({
    
    eta_state      <- cbind(0,fit1_dp_ranef_s$beta_RF)[,as.numeric(input$state)]
    eta_child      <- cbind(0,fit1_dp_ranef_s$beta_Fix[,1:2])[,as.numeric(input$child)]
    eta_area       <- cbind(0,fit1_dp_ranef_s$beta_Fix[,3])[,as.numeric(input$area)]
    eta_religion   <- cbind(0,fit1_dp_ranef_s$beta_Fix[,4:6])[,as.numeric(input$religion)]
    eta_education  <- cbind(0,fit1_dp_ranef_s$beta_Fix[,7:9])[,as.numeric(input$education)]
    eta_age        <- c(crossprod(B[input$age - 14,],t(fit1_dp_ranef_s$beta_spline)))
    
    prob1          <- plogis(eta_state + eta_child + eta_area + eta_religion + eta_education + eta_age)
    prob1  
  })
  
  prob2 <- reactive({
    
    eta_state      <- cbind(0,fit2_dp_ranef_s$beta_RF)[,as.numeric(input$state)]
    eta_child      <- cbind(0,fit2_dp_ranef_s$beta_Fix[,1:2])[,as.numeric(input$child)]
    eta_area       <- cbind(0,fit2_dp_ranef_s$beta_Fix[,3])[,as.numeric(input$area)]
    eta_religion   <- cbind(0,fit2_dp_ranef_s$beta_Fix[,4:6])[,as.numeric(input$religion)]
    eta_education  <- cbind(0,fit2_dp_ranef_s$beta_Fix[,7:9])[,as.numeric(input$education)]
    eta_age        <- c(crossprod(B[input$age - 14,],t(fit2_dp_ranef_s$beta_spline)))
    
    prob2          <- plogis(eta_state + eta_child + eta_area + eta_religion + eta_education + eta_age)

    prob2
  })
  
  prob3 <- reactive({
    
    eta_state      <- cbind(0,fit3_dp_ranef_s$beta_RF)[,as.numeric(input$state)]
    eta_child      <- cbind(0,fit3_dp_ranef_s$beta_Fix[,1:2])[,as.numeric(input$child)]
    eta_area       <- cbind(0,fit3_dp_ranef_s$beta_Fix[,3])[,as.numeric(input$area)]
    eta_religion   <- cbind(0,fit3_dp_ranef_s$beta_Fix[,4:6])[,as.numeric(input$religion)]
    eta_education  <- cbind(0,fit3_dp_ranef_s$beta_Fix[,7:9])[,as.numeric(input$education)]
    eta_age        <- c(crossprod(B[input$age - 14,],t(fit3_dp_ranef_s$beta_spline)))
    
    prob3 <- plogis(eta_state + eta_child + eta_area + eta_religion + eta_education + eta_age)
    
    prob3
  })
  
  
  output$UsageChoice <- renderPlot({
    ggplot(data=NULL,aes(x=prob1())) + geom_density(fill="blue",alpha=0.3) + xlab("Probability of using a contraceptive method") + ylab("Density")
  })
  
  output$ReversibleChoice <- renderPlot({
    ggplot(data=NULL,aes(x=prob2())) + geom_density(fill="blue",alpha=0.3) + xlab("Probability of using a reversible method") + ylab("Density")
  })
  
  output$MethodChoice <- renderPlot({
    ggplot(data=NULL,aes(x=prob3())) + geom_density(fill="blue",alpha=0.3) + xlab("Probability of using a modern method") + ylab("Density")
  })
  
  output$conditional <- renderTable({
    out <- data.frame(rho_1 = mean(prob1()),rho_2=mean(prob2()),rho_3 = mean(prob3()))
    colnames(out) <- c("Usage choice","Reversibility choice","Method choice")
    out
      })
  
  output$marginal <- renderTable({
    out <- data.frame(rho_1 = mean(1 - prob1()),
                      rho_2 = mean(prob1()*(1-prob2())),
                      rho_3 = mean(prob1()*prob2()*(1-prob3())),
                      rho_4 = mean(prob1()*prob2()*prob3()))
    colnames(out) <- c("No contraception","Sterilization","Natural methods","Modern methods")
    out
  })
  
})
