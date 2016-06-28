#
# Process data
#

process = function(DATA) {
hsA=colSums(DATA)
isA=rowSums(DATA)
colA=seq(1,length(hsA))
rowA=seq(0,length(hsA))

xA=sum(isA*rowA)
# Final size

N=sum(hsA*colA)
# Population size

list(hsA, isA, colA, rowA, xA, N)
}
