heter.f = function(loci.v) {
  
  loci.v = loci.v[loci.v!=0]  ## exclude the 0s
  if(length(loci.v)>0) {
  diff.v = diff(loci.v)[seq(1,length(diff(loci.v)),2)]
  sum.het = length(diff.v[diff.v!=0])
  sum.hom = length(diff.v[diff.v==0])
  sum.het/(length(diff.v)) } else {sum.het = NA}
  
}

#loci.v = zak.geno[1,6:(ncol(zak.geno))]
