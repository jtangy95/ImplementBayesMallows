sushi=read.table('/Users/changtaeyeong/Desktop/Mallows_rank_ model_Notes /sushi_Kamishima/sushi3a.5000.10.order', header=T, sep='\t')
head(sushi, 20)
str(sushi)
dim(sushi)

#install.packages('stringr')
library(stringr)

# Example of what 'str_match_all' function can do
test.input=('A4 pad notebook')
print(test.input)
Match<-str_match_all(test.input, "\\S+")
DF<-as.data.frame(Match)
names(DF)<-("Test")
print(DF)

# Use 'str_match_all' to convert string input(numbers with spaces) into numerical vector
input1=sushi[1, ]
Match<-str_match_all(input1, "\\d+")
v<-as.vector(Match)
output1=as.numeric(v[[1]][-(1:2),1])
print(output1)

N=5000  # Number of assessors
n=10  # Number of items

A<-matrix(0, ncol=N, nrow=n)
for(j in 1:N){
  input=sushi[j,]
  match=str_match_all(input, '\\d+')
  vec<-as.vector(match)
  output=as.numeric(vec[[1]][-(1:2), 1])
  A[,j]=output
}

print(A[, 1:10])
print(A[, (N-10):N])

# A does not exactly agree with R matrix in the paper.
# In R matrix...
# Row1 is item 0, Row2 is item 1, ... , and Row10 is item 9
# Use 'match' function to convert A into R matrix

output1
match(7, output1)

Rmat<-matrix(0, ncol=N, nrow=n)
for(j in 1:N){
  for(i in 1:n){
    Rmat[i,j]=match(i-1,A[,j])
  }
}
rownames(Rmat)<-paste('item', 0:9)
colnames(Rmat)<-paste('assessor', 1:N)

print(Rmat[,1:10])
print(Rmat[, (N-10):N])



