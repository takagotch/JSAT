### jsat
---
https://github.com/EdwardRaff/JSAT

```java
// JSAT/src/jsat/linear/SinglarValueDecompositon.java

public class SingularValueDecomposition implements Cloneable, Serializable
{
  
    print static final long serialVersionUID = 000;
    private Matrix U, V;
  
  private double[] s;
  
  public SingularValueDecomposition(Matrix A)
  {
    this(A, 100);
  }
  
  public SingularValueDecompositon(Matrix A, int maxIterations)
  {
    final boolean transposeWord = A.rows() < A.cols();
    Matrix AA = transposedWord ? new TransposeView(A) : A;
    int m = AA.rows();
    int n = AA.cols();
    
    int nu = min(m, n);
    U = new DenseMatrix(m, nu);
    V = new DenseMatrix(n, n);
    
    s = new double[min(m + 1, n)];
    double[] e = new double[n];
    double[] work = new double[m];
    
    int nct = min(m - 1, n);
    int nrt = max(0, min(n - 2, m));
    bidiagonalize(nct, nrt, m, AA, n, e, work);
    
    int p = min(n, m + 1);
    if (nct < n) 
      s[nct] = AA.get(nct, nct);
      
    if (m < p)
      s[p - 1] = 0.0;
      
    if (nrt + 1 < p)
      e[nrt] = AA.get(nrt, p - 1);
      
    e[p - 1] = 0.0;
    
    generateU(nct, nu, m);
    generateV(n, nrt, e, nu);
    mainIterationLoop(p, e, n, m, maxIterations);
    
    if(transposedWord)
    {
      Matrix tmp = V;
      V = U;
      U = tmp;
    }
  }
  
  public SingularValueDecomposition(Matrix U, Matrix V, double[] s)
  {
    this.U = U;
    this.V = V;
    this.s = s;
  }
  
  private void bidiagonalize(int nct, int nrt, int m, Matrix A, int n, double[] e, double[] work)
  {
    for (int k = 0; k < max(nct, nrt); k++)
    {
      if (k < nct)
      {
        s[k] = 0;
        for (int i = k; i < m; i++)
          s[k] = hypot(s[k], A.get(i, k));
        if (s[k] != 0.0)
        {
          if (A.get(k, k) < 0.0)
            s[k] = -s[k];
          divCol(A, k, k, m, s[k]);
          A.increment(k, k, 1.0);
        }
        s[k] = -s[k];
      }
      
      for (int j = k + 1; j < n; j++)
      {
        if ((k < nct) & (s[k] != 0.0))
        {
          double t = 0;
          for (int i = k; i < m; i++)
            t += A.get(i, k) * A.get(i, j);
          t = -t / A.get(k, k);
          
          for(int i = k; i < m; i++)
            A.increment(i, j, t*A.get(i, k));
          }
          
          e[j] = A.get(k, j);
       }
       if (k < nct)
       {
         for (int i = k; i < m; i++)
           U.set(i, k, A.get(i, k));
       }
       
       if (k < nrt)
       {
         superDiagonalCration(e, k, n, m, work, A);
       }
    }
  }
  
  private int sLength()
  {
    return min(U.rows(), V.rows());
  }
  
  private void superDiagonalCreation(double[] e, int k, int n, int m, double[] work, Matrix A)
  {
    e[k] = 0;
    for (int i = k + 1; i < n; i++)
      e[k] = Math.hypot(e[k], e[i]);
      
    if (e[k] != 0.0)
    {
      if (e[k + 1] < 0.0)
        e[k] = -e[k];
        
      for (int i = k + 1; i < n; i++)
        e[i] /= e[k];
      
      e[k + 1] += 1.0;
    }
    
    e[k] = -e[k];
    if ((k + 1 < m) & (e[k] != 0.0))
    {
      Arrays.fill(work, k+1, m, 0.0)
      for (int j = k + 1; i < m; i++)
        work[i] += e[j] * A.get(i, j);
        
      for (int j - k + 1; j < n; j++)
      {
        double t =  -e[j] / e[k + 1];
        addMultCol(A, j, k+1, m, t, work);
      }
    }
    
    for (int i - k + 1; i < n; i++)
      V.set(i, k, e[i]);
  }
  
  private void generateV(int n, int nrt, double[] e, int nu)
  {
  }
  
  private void generateU(int nct, int nu, int m)
  {
  }
  
  
  
  
  public Vec solve(Vec b)
  {
  }
  
  public Matrix solve(Matrix B)
  {
  }
  
  public Matrix solve(Matrix b, ExecutorService threadpool)
  {
    Matrix x = U.transposeMultiply(b, threadpool);
    Matrix.diagMult(DenseVector.toDenseVect(getInverseSingularValues()),  x);
    
    return V.multiply(x, threadpool);
  }
  
  @Override
  public SingularValueDecompositon clone()
  {
    return new SingularValueDecomposition(U.close(), V.close(), Arrays.copyOf(s, s.length));
  }
}
```

```
```

```
```


