{
InnerProd(b1,b2) =
  sum( i = 1, #b1, b1[i] * b2[i] );
}

{
  RandBasis(n,k) =
  matrix(n, n, x, y, random(2* k) - k);
}

{
GramSchmidt(B) =
  local(n,Bs,ip);
  n = #B;
  Bs = B;
  for( i = 1, n,
    for( j = 1, i-1,
      ip = InnerProd(B[,i], Bs[,j]) / InnerProd(Bs[,j] ,Bs[,j]);
      Bs[,i] = Bs[,i] - ip * Bs[,j]
    );
  );
return(Bs);
}

{
RoundOff(B, t) =
local(w,x);
  w = (B^-1) * t;
  x = round(w);
  return(B * x);
}

{
NearestPlane(B, Bs, t) =
local(n,v,f);
n = #B;
f = t;
v = 0 * t;
for(ii = 0, n-1,
      i = n -ii;
      y = InnerProd(f, Bs[,i]) / InnerProd(Bs[,i] , Bs[,i]);
      x = round(y);
      f = f - x * B[,i];
      v = v + x * B[,i];
  );
  return(v);
}
