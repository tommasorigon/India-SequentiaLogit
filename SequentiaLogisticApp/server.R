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
   
  output$UsageChoice <- renderPlot({
    eta_state      <- cbind(0,fit1_dp_ranef_s$beta_RF)[,as.numeric(input$state)]
    eta_child      <- cbind(0,fit1_dp_ranef_s$beta_Fix[,1:2])[,as.numeric(input$child)]
    eta_area       <- cbind(0,fit1_dp_ranef_s$beta_Fix[,3])[,as.numeric(input$area)]
    eta_religion   <- cbind(0,fit1_dp_ranef_s$beta_Fix[,4:6])[,as.numeric(input$religion)]
    eta_education  <- cbind(0,fit1_dp_ranef_s$beta_Fix[,7:9])[,as.numeric(input$education)]
    eta_age        <- c(crossprod(B[input$age - 14,],t(fit1_dp_ranef_s$beta_spline)))
    
    prob <- plogis(eta_state + eta_child + eta_area + eta_religion + eta_education + eta_age)
    # Final plot
    ggplot(data=NULL,aes(x=prob)) + geom_density(fill="blue",alpha=0.3) + xlab("Probability of using a contraceptive method") + ylab("Density")
    
  })
  
  output$ReversibleChoice <- renderPlot({
    
    eta_state      <- cbind(0,fit2_dp_ranef_s$beta_RF)[,as.numeric(input$state)]
    eta_child      <- cbind(0,fit2_dp_ranef_s$beta_Fix[,1:2])[,as.numeric(input$child)]
    eta_area       <- cbind(0,fit2_dp_ranef_s$beta_Fix[,3])[,as.numeric(input$area)]
    eta_religion   <- cbind(0,fit2_dp_ranef_s$beta_Fix[,4:6])[,as.numeric(input$religion)]
    eta_education  <- cbind(0,fit2_dp_ranef_s$beta_Fix[,7:9])[,as.numeric(input$education)]
    eta_age        <- c(crossprod(B[input$age - 14,],t(fit2_dp_ranef_s$beta_spline)))
    
    prob <- plogis(eta_state + eta_child + eta_area + eta_religion + eta_education + eta_age)
    # Final plot
    ggplot(data=NULL,aes(x=prob)) + geom_density(fill="blue",alpha=0.3) + xlab("Probability of using a reversible method") + ylab("Density")
    
  })
  
  output$MethodChoice <- renderPlot({
    
    eta_state      <- cbind(0,fit3_dp_ranef_s$beta_RF)[,as.numeric(input$state)]
    eta_child      <- cbind(0,fit3_dp_ranef_s$beta_Fix[,1:2])[,as.numeric(input$child)]
    eta_area       <- cbind(0,fit3_dp_ranef_s$beta_Fix[,3])[,as.numeric(input$area)]
    eta_religion   <- cbind(0,fit3_dp_ranef_s$beta_Fix[,4:6])[,as.numeric(input$religion)]
    eta_education  <- cbind(0,fit3_dp_ranef_s$beta_Fix[,7:9])[,as.numeric(input$education)]
    eta_age        <- c(crossprod(B[input$age - 14,],t(fit3_dp_ranef_s$beta_spline)))
    
    prob <- plogis(eta_state + eta_child + eta_area + eta_religion + eta_education + eta_age)
    # Final plot
    ggplot(data=NULL,aes(x=prob)) + geom_density(fill="blue",alpha=0.3) + xlab("Probability of using a traditional method") + ylab("Density")
    
  })
  
})
