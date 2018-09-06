                                                   

                                                            #GGRK!  HK! LIMTG#

#####R revision!#####

#Working directory
getwd()
setwd("E:/R")

#Packages
install.packages('PackageName') #install package. Note use of " "
library(PackageName) # Open package for use...note " " not required now
detach(package=PackageName) # Close package so that it does not interfere with other packages

install.packages("installr"); 
library(installr) # install+load installr


#Most basic unit: Vector! (Like Matrices in MATLAB)

remove(a)   #removing variables
rm(list=ls( )) #clear all      #ls() is like whos() in MATLAB
#Note: Just typing the function name without parenthesis gives the complete code of the function in R
dev.off() # removes all graphics

?a( ) or help('a') #for help
?"+"  #needs to be used when asking for help with symbols

#One less discussed class is Integers
a<-c(3L,5L) #adding L makes it integer class
class(a)
#integer class consists of no decimals and hence easier to process than numeric class...particularly for large data sets

seq(1,10,2) #generates vectror from 1 to 10 in steps of 2...default third argument is the step size
seq(1,10,length.out=100) # you can specify an argument....length.out here determines the vector length

Length(a)  Class(a)

data("dataset") #calls the dataset as variable   and str(dataset) tells you about the variables and their numbers. Names( ) can be used to reveal just the column names. Note that once the variable comes into the workspace, no need for inverted commas!
dataframe$ColumnName as well as dataframe[["columnName"]] works! 
Dataframe["columnName"] also works, but the class is a data frame, not a vector!
  
Head( ) and tail( ) can be used to view parts of data!

$ sign #is the 'accessor' used to access data through column heads

#Factor variables store categorical data!
levels(factorVariable) #gives you different categories

identical() # use to check if two vectors are equal in all ways


table( ) #is typically used to get the frequencies values of a vector

#Writing words inside vectors with an equal to sign makes it a column name!
a<-c(Italy=100,Germany=300)   #or
names()#can be used to assign or read column names

#Subsetting! : access specific parts of Vectors
a[1]; a[c(1,3)]
a["ColumName"]; a[c("ColName1","ColName2")] #Note that indexing directly accepts column names.


#Coercion a unique feature of R
#Many times, it tries to imagine what you meant to code and try to process it that way!!!
a<-c(1,"Hi",2) #must ideally show an error since a vector must have similar class entries and here, we have a chaacter entry with numerics!
# R coerces even 1,2 to characters in this case and assigns this vector's class as character!!! Try out---
a
class(a)
#Coercion can lead to missing values or NA (as R designates it)
a<-c("1","hi","2")
as.numeric(a)  #leads to a NA in the second place since it could not coerce "hi" into a numeral! 


sort() # arranges data in ascending order (default)
order() #does the same with indexes..gets the index numbers of values in ac=scending order...useful if you think about it...
which.max() #gives the index of the highest value
rank() #ranks each value based on its magnitude in the vector


#DataFrame!  Can be used to convert vectors into columns...note that column names can be specified!
temp <- c(35, 88, 42, 84, 81, 30)
city <- c("Beijing", "Lagos", "Paris", "Rio de Janeiro", "San Juan", "Toronto")
city_temps <- data.frame(name = city, temperature = temp)
#Note however that charcters are converted into factors by default: class(city_temps$name)...to avoid this..
city_temps <- data.frame(name = city, temperature = temp, stringsAsFactors=FALSE)


#NA is what appears for missing values in R
#They must be identified and removed before processing data since descriptive o/p will be NA till all NAs are removed!
is.na(a) # logical vector identifies indexes with NA
#Use the negation operator !  to remove NA
a<-a[!is.na(a)]



#Using logical vectors to index things is a beautiful thing in R
a<-c(1,2,3,5,9,10,14,17,19,34)
index<-a>5
b<-a[index] #gives me only those >5
#Extending on this, suppose I want >5 and is odd....create the two logical vectos and index them to a!
G5<-a>5
odd<-a/2!=round(a/2)#dividing by 2 and identifying odd numbers by decimal answers!
b<-a[G5&odd]


which()#can be used to get the indices of a logical vector that are true...this is not strictly needed since a logical vector can directly be used for indexing
#  ie takes a logical vector and tells you which indices are true
#but using which() reduces the processing load
match() #can be used when there are many items in a vector to compare with the second vector
a<-c(1,2,3,5,9,10,14,17,19,34)
b<-c(1,14,34)
match(b,a) #gives you indices of a which matched with b!
#Suppose what you are matching is not there in the second vector?
c<-c(5,9,30) #30 is not there in a 
match(c,a) #you will get Na for 30 indicating such number doesnt exist in a
%in% # can also be used to check if value sin b are present in a
b %in% a  #gives a logical vector for each position!


#Using the dplyr package!
 # It is a very powerful package for datatable manipulation. More efficient and intuitive than
 #base R functions.
mutate(dataframe,rate=value/pop)# adds another column-'rate' to the dataframe using the existing columns value and pop in the dataframe. Note that dataframe need not be referred to again and again!
filter(datframe,rate>0.5)#applies the logical and gives a filtered dataframe...using base functions, we needed to have a index vector of logicals and subset it to the dataframe...this is shorter
select(dataframe,value,pop)#gives only those columns in the new table
#Another strong feature is Piping! : %>%
 #Mention dataframe once and carry out operations as if you are writing about it! Very inntuitive!
dataframe %>% select(C1name,C2name,C3name) %>% filter(C3name>0.7)
my_states<-murders %>% mutate(rate=(total/population)*100000,rank=rank(-rate))%>% filter(region %in% c("Northeast","West") & rate<1) %>% select(state,rate, rank)

#Exploratory data analysis/Plots- the strength of R!!!
x<-c(1,2,3); y=c(5,6,7)
plot(log(x),log(y)) #Put what you want on the axis labels inside plot....define them before you put it into plots for that convenience
hist(x) #histogram
boxplot(y~x,data=murders) # boxplot...when you use tilda, y axis or dep variable comes first and X axis/IV comes next...note that dataset needs to be mentioned!
lines() # is a sort of base plotting building over previous plot like hold on
a<-c(10,22,30); b<-c(14,15,16); c<-c(12,14,16);d<-c(3,15,7);
plot(a,b)
lines(c,d)




                   #Conditionals!!!...building blocks of R programming!

   #Condition is housed in small brackets, action come next in flower brackets...note that the initial flower bracket starts with the conditional itself...and ends with a lone flower bracket!
a<-c(1,2,3)
if (a!=0){
  print(1/a)
} else {
  print("No reciprocal for 0")
}
#Put the most frequent possibility first so that it need not go to else in the first place! Speeds things up!


#Another related function is ifelse()
#It takes 3 arguments- a logical, and two answers...if the logical is true, first ans is printed and if it is false, the second one if printed.
a<-0
ifelse(a>0,1/a,NA)

a<-c(1,2,4,NA,8,10,NA,89)
No_Nas<-ifelse(is.na(a),0,a) #Replacing NAs with Zero in vector a
No_Nas
#To remove NAs completely, you could just use logical vector indexing!
a[!is.na(a)]

any() #takes a logical vector and returns true if there is atleast one true!
all() #takes a logical vector and returns true if and only if all are true!


          ###Functions!!!###
a<-c(1,2,3,4,5)
#Making a function...
avg<-function(x){
  s<-sum(x)
  l<-length(x)
  s/l
}
#Using it...
avg(a)
#Functions can have more than one variable...function(x,y,z)
avg<-function(x,arithmetic=TRUE){
  n<-length(x)
  ifelse(arithmetic,sum(x)/n,prod(x)^(1/n))
}
#Here, if arithmetic is false, geometric mean is the output!

   # For-loops  (They use Empty Vector-Indexing)
#Say, we need to know the sum of objects till a particular number
compute_sum<-function(n){
  x<-1:n
  sum(x)
}
compute_sum(5)
#What if we need to do this for all numbers from 1 till 25?...Same thing is being done for a range of values of n
Compute_range_sum<-function(x){
  #create empty vector
  m<-vector(length=x)
  #Apply for-loop and index empty vector to n
  for (n in 1:x){
    m[n]<-sum(1:n)
  }
  m #Note that m is coming after the for-loop, not inside!
}
Compute_range_sum(7)



# However, we have stronger and more efficient functions than loops in R. They are...
#The apply family!
apply()
sapply()
tapply()
mapply()
#Other important funstions include-
split(); cut(); reduce(); quantile(); unique(); identical()



                                ###  DATA VISUALIZATION  ###
table() # can be used to get the frequency information of each point of data...like a histogram basically.
#Cumulative distribution...F(a)=Pr(x<=a)...ie the proportion/probability of values lesser than or equal to a value.
#Histogram is a better descriptor of distribution of data! However, we lose info about data points falling in the same bin.
#Smooth density curves just plot the envelope of the histogram...also, Y axis has density/frequency rather than counts.
 #Area under the curve is 1...this is aesthetically more pleasing than hist and comparisons between distributions are easier.(standardized)
 #Comparable to a normal distribution...the two can be overlapped to see if it fits into normal distribution.

#Standard units: (x-mean)/SD  The standard units immedietely tell you how far away a value is from the mean...also note that original units do not matter.
scale(x) # can be used to convert a value into standard units
mean(abs(z)<2) # is around 0.95 if the distribution is normal...so convert data into std units n see the area within 2 SD...similar to 95% means distribution is Gaussian
#or better, just overlap on normal distribution and see!

#Once we understand that a distru=ibution is nearly normal, just two measures- Mean and SD can be used to model the distribution!
norm() #series can be used...like pnorm(): CDF,qnorm : inverse of CDF, dnorm(): density distribution and rnorm(): generate numbers from this set
#Note though that the approximation may not be very accurate at the tails....the middle values can be more safely predicted.
pnorm(72, mean=70, sd=5)-pnorm(69, mean=70, sd=4) # gives us the proportion of values between 69 and 72 in a normal distribution centred around 70. Note though that this is around centre of distribution and can be quite safely predicted.
pnorm(81, mean=70, sd=5)-pnorm(79, mean=70, sd=4) #are the tail values and hence are more difficult to model correctly. Normal distribution tends to underestimate tail values...ie fatter tails are not predicted well...extreme values are lesser than what they actually are. 


#Quantiles are the other counterpart of cumulative probabilities
#For a given p, the counterpart on the x axis of the distribution is the Quantile (q)
#For eg, for a p=0.5, q can be 69...if it is a normal distribution, the quantiles will overlap
p<-seq(0.05, 0.95, 0.05 )
observed_quantiles<-quantile(x,p) # where X is the vector for which quantiles are being sought
theoretical_qunatiles<-qnorm(p,mean=mean(x),sd=sd(x))
plot(observed_quantiles,theoretical_qunatiles)
abline(0,1)# draws a straight line with a intercept of zero and slope of 1...ie a 45 deg linear line
#Converting ito z units makes things easier with no need for mean, SD in the theoretical quantiles
z<-scale(x) 
observed_quantiles<-quantile(z,p)
theoretical_quantiles<-qnorm(p) #....and then plot
#Percentiles are special case of quantiles....for p<-seq(0.01,0.99,0.01) ...ie 1% to 99%....special case are  quartiles: 0.25, 0.5, 0.75


                             ###Plotting with GGPlot2!!!### (Works with Data Tables only!)

#GGPlot or Grammar of Graphics + Plot Has 5 components to it- 
#Data, Geometry(type of plot), Mapping (what goes where), Scale and Labelling

#Data component
ggplot(data=murders)
#or pipe in
murders %>% ggplot()
#You can assign an object so that plot immediately doesn't happen
p<-ggplot(data=murders) # Print it when needed with print(p) or just p

#ggplot hence has layers to it concatenated by +
Data %>% ggplot()+layer2+l3+....   #or
ggplot(data=Data)+l2+.....

#Adding Geometry component
murders %>% ggplot()+geom_point(aes(x=pop,y=total)) #aes() is the function to used to locate the variables mentioned, and as such can be found in every layer! It connects data table with the graph!
        #geom_point() refers to a scatter plot
#Adding labels to each point on scatter plot
murders %>% ggplot()+geom_point(aes(x=pop,y=total)) +geom_text(aes(x=pop,y=total,label=abb)) # geom_text is being given the info about each data point position and which variable in the data table to take labelling from.

#or you can put aesthetic mapping first...and geometry later...this way, global mapping takes over and we need not type things again and again
p<-murders %>% ggplot(aes(x=pop,y=total, label=abb))+geom_point(size=3)+geom_text(nudge_x=3) #nudge_x here spaces the labels such that they don't overlap on each other

#making the scale log
p<-p+scale_x_log10()+scale_y_log10() #might have to adjust nudge_x value now that axis has become log

#add labels
p+xlab(".....")+ylab(".....")+ggtitle(".......")

# Changing color of points...
geom_point(color="blue",size=3) #makes the color blue
geom_point(aes(color=region),size=3) # assigns color based on region...note that aes() is used when a variable is mentioned..this also adds a legend by default.

# More changes can be made regarding other finer aspects like Legend text, adding abline, standard themes etc
#For abline, typically, intercept=0, slope=1...often, you out the log10(average) as intercept 
r<-murders %>% summarise(rate=sum(total)/sum(population)*10^6) %>% .$rate # .$rate converts r into a number rather than being a dataframe (somewhat like a virtual column?)
geom_abline(intercept=log10(r))


#Say we want to plot a histogram

p<-heights %>% filter(sex="Male") %>% ggplot(aes(x=height)) # Only one axis required for histogram
p1<-p+goem_histogram(binwidth=1,fill="blue",color="black")
p2<-p+goem_histogram(binwidth=2,fill="blue",color="black")
p3<-p+goem_histogram(binwidth=3,fill="blue",color="black")
library(gridExtra)
grid.arrange(p1,p2,p3, ncol=3) #Plotting many plots in one



#Libraries important for plotting!
library(dplyr)
library(ggplot2)
library(dslabs)
library(tidyverse)

data(murders)
murders %>% ggplot(aes(population, total, label = abb, color = region))+geom_point() #note that even though label argument is given, it is not executed since there is no geom_label() command
murders %>% ggplot(aes(population, total, label = abb, color = region))+geom_label()+geom_point() # here, the command is executed

data(heights)
#Two plots within the same figure by using group argument in ggplot()
heights %>% ggplot(aes(height, group = sex)) + geom_density()
#You can change fill which automatically assumes many groups
heights %>% ggplot(aes(height, fill = sex)) + geom_density()
#to avoid overlapping of density plots, use alpha bending...ie use the alpha argument in geom_density() to produce transparency
heights %>% ggplot(aes(height, fill = sex)) + geom_density(alpha=0.2) 


# Summarise() function
library(tidyverse)
library(dslabs)
data(heights)
# Summarise() can be used to directly access variables in the data frame. Also, it forms a dataframe of the summary variables in the argument which can be easily accessed with the $ sign
s<-heights %>% filter(sex=="Male") %>% summarise(average=mean(height), SD=sd(height)) 
s$average
# similarly
s<-heights %>% filter(sex=="Male") %>% summarise(average=mean(height), SD=sd(height),minimum=min(height),Maximum=max(height)) # can keep on adding
# a caveat though....summarise() gives each variable a single value only...those with multiple values give an error....eg: summarise(range=quantile(vector,q=seq(0.1,0.9,0.1)))


# dot placeholder
#Note that most of the output in dplyr package are in dataframes....if we need to use the output as an input to the function which needs a numeric, can be used.
data(murders)
US_murder_rate<- murders %>% summarize(rate=sum(total)/sum(population)*100000)
# US_murder_rate is a datframe...to make it a numeric, use the dot placeholder
b<-murders %>% summarize(rate=sum(total)/sum(population)*100000) %>% .$rate
class(b)


#group_by() function
#Often used to separate out data into different parts before summarise()
a<-heights %>% group_by(sex) # converts into different groups internally...even the class actually changes into a group dataframe
# So, basically a preceeding opearator to other functions...like 'split data' in SPSS
class(a) # summarise() behaves differently after applying group_by()
heights %>% group_by(sex) %>% summarise(mean=mean(height),SD=sd(height))

# arrange() function is very useful in dplyr package to arrange data
murders %>% arrange(population) %>% head() # arrange data by asc order of pop
murders %>% arrange(desc(population)) %>% head()  # descending order
# Can use two arguments too....ordering within the clustered data...
murders %>% arrange(region,population) %>% head()
#top_n() function to see the top 10 or so (instead of head())
murders %>% arrange(region,population) %>% top_n(10)


                                #Gapminder dataset analysis!!!

library(dslabs)
data(gapminder)
head(gapminder)

#Looking at infant mortality in two countries in the year 2015
gapminder %>% filter(year==2015 & country %in% c("Sri Lanka","Turkey")) %>% select(country, infant_mortality) #Note that filter() takes a logical input

ds_theme_set() # particular theme about background etc
gapminder %>% filter(year==1962) %>% ggplot(aes(fertility,life_expectancy,col=continent)) + geom_point()
#One sees that Africa and Asia have high fertility and low life expectancy

#Does it hold in 2012?  FACETING-
#facet_grid() can be used to split the figure into rows and columns
gapminder %>% filter(year %in% c(1962,2012)) %>% ggplot(aes(fertility,life_expectancy,col=continent))+
geom_point()+facet_grid(continent~year)
#Want stratification just by year? Use dot place holder
gapminder %>% filter(year %in% c(1962,2012)) %>% ggplot(aes(fertility,life_expectancy,col=continent))+
geom_point()+facet_grid(.~year)
#Looking just at Europe and Asia...and across years: facet_wrap is useful when there are many plots
gapminder %>% filter(year %in% c(1962,1980, 1990,2000,2012) & continent %in% c("Europe","Asia")) %>% ggplot(aes(fertility,life_expectancy,col=continent))+
geom_point()+facet_wrap(~year) #Not dot place-holder necessary for facet_wrap() since it is built to handle many variables (not exactly two variables as in dacet_grid())

#Comparing two countries-
countries<-c("India","Germany")
gapminder %>% filter(country %in% countries) %>% ggplot(aes(year,fertility,group=country))+geom_line() #group argument is required for two separate lines!
#Using col argument would automatically group and add a legend!
gapminder %>% filter(country %in% countries) %>% ggplot(aes(year,fertility,col=country))+geom_line()
#Rather than using legends, labelling directly on plot is better!
#For this, first create a datframe with the labels and their x and y positions
Label<-data.frame(country=countries,x=c(1975,1965),y=c(4.5,2)) #x axis is in years, y axis is fertility
#Now add geom_label()
gapminder %>% filter(country %in% countries) %>% ggplot(aes(year,fertility,col=country))+geom_line()+
geom_label(data=Label,aes(x,y,label=country),size=5)
#Remove legend with theme()
gapminder %>% filter(country %in% countries) %>% ggplot(aes(year,fertility,col=country))+geom_line()+
geom_label(data=Label,aes(x,y,label=country),size=5)+theme(legend.position="none")

#Transformations
#You can transform axis or transform the data!
#Transforming axis for example-
gapminder %>% filter(country %in% countries) %>% ggplot(aes(year,fertility,col=country))+geom_line()+
geom_label(data=Label,aes(x,y,label=country),size=5)+scale_x_continuous(trans="log2") #scale_x_continous is used to change axes
#...or transform data at the ggplot() level itself as in log2(year)

    #Stratify and boxplot
#include removal of NAs in filter()
p<-gapminder%>%mutate(dollars_per_day=gdp/population/365)%>%
filter(year==2000 & !is.na(gdp))%>%ggplot(aes(region,dollars_per_day))+geom_boxplot()
#theme() can be used to make changes in figures like- axis text
p+theme(axis.text.x=element_text(angle=90,hjust=1)) #text made vertical, with a justification of 1
#we can reorder text according to say -median of dollars_per_day instead of alphabetical
p<-gapminder%>%mutate(dollars_per_day=gdp/population/365)%>%filter(year==2000 & !is.na(gdp))%>%
mutate(region=reorder(region,dollars_per_day,FUN=median))%>%
ggplot(aes(region,dollars_per_day,fill=continent))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1))
#reorder() basically reorders a factor variable in accordance to an operation/function of another variable
#Note: filter() #often used to remove NAs from the groups that are being formed
p+scale_y_continuous(trans="log2")+geom_point(show.legend=FALSE) #a log transform of y axis + showing all data!

intersect() #can be used to take only the ones which are common to two vectors
fill=factor(varible) # within ggplot(), this can be used to split data and plot...factor() can be very useful like this
case_when() #very useful when choosing to address case by case basis. Example-
gapminder %>% mutate(group=case_when(.$region %in% "west"~"The West"+.$region %in% c("Asia","Southeast Asia")~"Asia")) # etc

#Note- fill,labels etc must be used inside aes()...they need access to variables you see.










































