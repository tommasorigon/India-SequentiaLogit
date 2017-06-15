## ui.R ##
library(shinydashboard)

header <- dashboardHeader(
  title = "Sequential probabilities"
)

body <- dashboardBody(
  fluidRow(
    column(width = 7,
           box(width = NULL,
               h4("Usage Choice"),p("In the graph below is represented the posterior density of the probability of using a contraceptive method."),
               plotOutput("UsageChoice")
           ),
           box(width = NULL,
               h4("Reversible Choice"),p("In the graph below is represented the posterior density of the probability of using a reversible method, provided that a contraceptive method was used."),
               plotOutput("ReversibleChoice")
           ),
           box(width = NULL,
               h4("Method Choice"),p("In the graph below is represented the posterior density of the probability of using a traditional method, provided that a reversible contraceptive method was used."),
               plotOutput("MethodChoice")
           )
    ),
    column(width = 5,
           
           box(width = NULL, status = "warning",
               h4("Sequential probabilities"),
               p("The values below represent the posterior estimate of the conditional probabilities associated to the sequential process."),
               tableOutput("conditional"),
               h4("Marginal probabilities"),
               p("The values below represent the posterior estimate of the probabilities of the multinomial sampling. Their summation is equal to one."),
               tableOutput("marginal")
           ),
           box(width = NULL, status = "warning",
               selectInput("state", label = h4("State"), 
                           choices = list("Andhra Pradesh"=2,"Arunachal Pradesh"=3, "Assam"=4,    
                                          "Bihar" =5, "Chandigarh" =6, "Chhattisgarh" =7,
                                          "Dadra+Nagar Haveli"=8,  "Daman & Diu" =9,      
                                          "NCT of Delhi"=10, "Goa" =11, "Gujarat" =12,           
                                          "Haryana"=13, "Himachal Pradesh"=14, "Jammu & Kashmir"=15,
                                          "Jharkhand"=16, "Karnataka"=17, "Kerala"=18,           
                                          "Madhya Pradesh"=19, "Maharashtra"=20, "Manipur"=21,         
                                          "Meghalaya"=22, "Mizoram"=23, "Nagaland"=24,         
                                          "Orissa"=25, "Pondicherry"=26, "Punjab"=27,            
                                          "Rajasthan" =28,"Sikkim"=29, "Tamil Nadu"=30,        
                                          "Tripura"=31, "Uttarakhand"=32,"Uttar Pradesh" = 1,
                                          "West Bengal"=33), selected = 1)
           ),
           box(width = NULL, status = "warning",
               sliderInput("age", label = h4("Age"), min = 15, 
                           max = 49, value = 30)
               
           ),
           box(width = NULL, status = "warning",
               radioButtons("child", label = h4("Child"),
                            choices = list("None" = 2, "One" = 3, "More than one" = 1), 
                            selected = 1),
               radioButtons("area", label = h4("Area"),
                            choices = list("rural" = 1, "urban" = 2), 
                            selected = 1),
               radioButtons("religion", label = h4("Religion"),
                            choices = list("Hindu" = 1, "Muslim" = 2, "Christian" = 3,"Other" = 4), 
                            selected = 1),
               radioButtons("education", label = h4("Education"),
                            choices = list("None" = 1, "Low" = 2, "Intermediate" = 3,"High" = 4), 
                            selected = 1)
           )
    )
  )
)

dashboardPage(
  header,
  dashboardSidebar(disable = TRUE),
  body
)