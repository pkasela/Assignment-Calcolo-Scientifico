%
%   [bp, bm] = bern(x)
%   
%    This function computes the inverse of the Bernoulli function 
%
%    bp = B(x)  = x / (exp(x)-1)
%    bm = B(-x) = x + B(x)
%
function [bp, bn] = bern (x)

  xlim = 1e-2;
  ax   = abs (x);


  % Calcola la funz. di Bernoulli per x=0

  if (ax == 0)
    bp=1.;
    bn=1.;
    return
  end;


  % Calcola la funz. di Bernoulli per valori
  % asintotici dell'argomento

  if (ax > 80),
    if (x >0),
      bp=0.;
      bn=x;
      return
    else
      bp=-x;
      bn=0.;
      return
    end;
  end;


  % Calcola la funz. di Bernoulli per valori
  % intermedi dell'argomento

  if (ax > xlim),
    bp=x./(exp(x)-1);
    bn=x+bp;
    return
  else

    
    % Calcola la funz. di Bernoulli per valori
    % piccoli dell'argomento mediante sviluppo
    % di Taylor troncato dell'esponenziale
    
    ii=1;
    fp=1.;
    fn=1.;
    df=1.;
    segno=1.;
    while (abs(df) > eps),
      ii=ii+1;
      segno=-segno;
      df=df*x/ii;
      fp=fp+df;
      fn=fn+segno*df;
      bp=1./fp;
      bn=1./fn;
    end;
    return
  end
