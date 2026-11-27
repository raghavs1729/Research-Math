# Direct numerical test of the congruences of arXiv:2603.08510 (closing remarks and Conjectures 1, 2).
import flint, sys
def kron(a,p):
    a%=p
    return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
def pbar_mod(ell, D):
    c=[0]*D; n=0
    while n*n<D: c[n*n]=1 if n==0 else (2*(-1)**n)%ell; n+=1
    return flint.nmod_poly(c,ell).inverse_series_trunc(D)   # 1/phi(-q) = sum pbar(n) q^n
D=int(sys.argv[1])
for (ell, mult, cond, desc) in [
    (17, 17*121, lambda n: n%8==3 and kron(n,11)==-1, 'pbar(17*11^2 n) mod 17, n=3 (8), (n/11)=-1'),
    (23, 23*169, lambda n: n%8==5 and kron(n,13)==1, 'pbar(23*13^2 n) mod 23, n=5 (8), (n/13)=+1'),
    (7, 16, lambda n: n%8==3 and kron(n,7)==1, 'pbar(16 n) mod 7, n=3 (8), (n/7)=+1'),
    (19, 4*289, lambda n: n%8==3 and kron(n,17)==-1 and kron(n,19)==1, 'pbar(4*17^2 n) mod 19, n=3 (8), (n/17)=-1, (n/19)=+1')]:
    pb=pbar_mod(ell,D)
    ns=[n for n in range((D-1)//mult+1) if cond(n)]
    bad=[n for n in ns if int(pb[mult*n])!=0]
    print(f'{desc}: tested {len(ns)} values of n (n <= {ns[-1] if ns else None}); nonzero: {bad[:10]}')
