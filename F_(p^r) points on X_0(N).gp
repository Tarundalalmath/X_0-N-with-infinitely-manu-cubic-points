\\ tr.gp 
\\ William A. Stein
\\ Created:   October 1998
\\ Modified: 14 January 2001
\\
\\ DESCRIPTION: Hijikata trace formula for tr T_n 
\\ on S_k(Gamma_0(N)), any n>=1, except if n|N, 
\\ then n must be prime.

{h(d)=
   \\ qfbclassno(-236) breaks, so use qfbclassno(-236,1), 
   \\ i.e., use Euler products instead!
   qfbclassno(d,1);
}

{w(d)=
   if(d==-4,4,if(d==-3,6,2))/2;
}

{tof(a, f=0,t=1,i=1,b=0)=
   f=factor(a);
   b=prod(i=1,matsize(f)[1],
          t*=f[i,1]^((f[i,2]-f[i,2]%2)/2)
     );
   if((a/t^2) % 4 == 1, t, t/2);
}

{mu(N, f)=
   f = factor(N); 
   N*prod(i=1,matsize(f)[1],(1+1/f[i,1]));
}

{sig(n, N, f=0)=
   if(N%n==0, 
      n, 
   \\ else
      f=factor(n);  
      prod(i=1,matsize(f)[1], 
         (1-f[i,1]^(f[i,2]+1))/(1-f[i,1])));
}

{quadpoints(s,n,p,v, len=0,pv=0,x=0,vec=0)= 
   len=0;
   pv=p^v;
   vec=vector(pv,x,0);
   for(x=0,pv-1,        \\ very simplistic!
      if((x^2-s*x+n)%pv==0,
         len++;
         vec[len]=x
      )
   ); 
   vector(len,x,vec[x]);
}

{A(s,f,n,p,v, sr=0, rho=0, cnt=0)=
   rho=valuation(f,p);
   sr=eval(Set(quadpoints(s,n,p,v+2*rho)%(p^(v+rho))));  
   cnt=0;
   for(i=1,length(sr),
      if((2*sr[i]-s)%(p^rho)==0 && if(n==p,sr[i]%p!=0,1), 
         cnt++
      )
   );
   cnt;
}

{B(s,f,n,p,v)=
   rho = valuation(f,p);
   length(Set(quadpoints(s,n,p,v+2*rho+1) % p^(v+rho)));
}

{cp(s,f,n,p,v)=
   A(s,f,n,p,v) +
      if(((s^2-4*n)/f^2)%p==0,
         B(s,f,n,p,v)
      );
}

{c(s,f,n,N, fac)=
   fac=factor(N);
   prod(i=1,matsize(fac)[1],
      cp(s,f,n,fac[i,1],fac[i,2])
   );
}

{type_p(n,k,N, s=0)=
   if(issquare(n), 
      s=floor(sqrt(4*n));
      1/4*(s/2)*n^(k/2-1)*(c(s,1,n,N) + (-1)^k*c(-s,1,n,N)),0
   );
}

{absxy(s,n,k,N,f=0,x=0,y=0,t=0)=
   t = floor(sqrt(s^2-4*n));
   x = (s-t)/2;
   y = (s+t)/2;
   (min(abs(x),abs(y))^(k-1)/abs(x-y)) 
       * sign(x)^k 
       * sumdiv(t,f,1/2*eulerphi(t/f)*c(s,f,n,N));
}

{type_h(n,k,N, s=0,sm=0)=
   for(s=ceil(2*sqrt(n))+issquare(n),4*n,
      if(issquare(s^2-4*n),
         sm += (absxy(s,n,k,N)+absxy(-s,n,k,N))
      )
   ); 
   sm;
}

{xy(s,n,k,N, f=0, x=0, y=0)=
   x = (s+sqrt(s^2-4*n))/2;
   y = (s-sqrt(s^2-4*n))/2;
   1/2 * (x^(k-1)-y^(k-1))/(x-y)
       * sumdiv(tof(s^2-4*n),f,
            h((s^2-4*n)/f^2)/w((s^2-4*n)/f^2) * c(s,f,n,N)
         ); 
}

{type_e(n,k,N,r=0)=
   r=floor(2*sqrt(n))-issquare(n);
   sum(s=-r,r,xy(s,n,k,N));
}

{sum_s(n,k,N)=
   type_p(n,k,N)+type_h(n,k,N)+type_e(n,k,N);
}

{tr(n,k=12,N=1, t=0)=
   t=(if(k%2,0,-sum_s(n,k,N)+if(k==2,sig(n,N))
         +if(issquare(n),(k-1)*mu(N)/12*n^(k/2-1))));
   if(abs(t-round(t))>.01,print("Precision loss -- Error! noninteger trace! ",t),round(t)); 
}

print("tr(n,k,N) = tr(T_n on S_k(Gamma_0(N))). (If n | N then n must be prime.)");
X(p,r,N)= 1+p^r-tr(p^r,2,N)+p*tr(p^(r-2),2,N);
print("X(p,r,N)= F_(p^r) points on X_0(N).");
