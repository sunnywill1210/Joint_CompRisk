### Data cleaning for BMT data
#bmt <- read.csv("http://cedric.cnam.fr/~latoucha/R/data.bmt.csv",stringsAsFactors = FALSE)
bmt <- read.csv(" .\\bmt_data.csv",stringsAsFactors = FALSE)

bmt$time <- bmt$SURVM
bmt$ev <- bmt$CRELAPSE

# 24 months followup 
bmt$time <- pmin(24, bmt$time)
bmt$ev <- ifelse(bmt$time == 24, 0, bmt$ev)

time=bmt$time
type=bmt$ev
adv =bmt$ADV
mini =bmt$MINI

data<-data.frame(time,type, adv, mini)
## time: event time
## ev: event type: 0 censored observation; 1 relapse; 2 TRM
## mini: Treatment arm: 0 MAC; 1 RIC
## adv:  Status at transplantation: 1 advanced; 0 Other



