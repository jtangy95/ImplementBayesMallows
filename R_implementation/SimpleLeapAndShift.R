n=8; L=3
r=sample(1:n, n)
print(r)
u=sample(1:n, 1)
print(u)
r.star=r
r.star[u]=sample(setdiff(max(r[u]-L,1):min(r[u]+L,n),r[u]),1)
print(r.star)
r.prime=r
delta=r.star[u]-r[u]
if(delta>0){
  r.prime[r[u]<r & r<=r.star[u]]=r.prime[r[u]<r & r<=r.star[u]]-1
} else{
  r.prime[r[u]>r & r.star[u]<=r ]=r.prime[r[u]>r & r.star[u]<=r ]+1
}
r.prime
r.prime[u]=r.star[u]
print(r) ; print(r.prime)

rm(list=ls())
