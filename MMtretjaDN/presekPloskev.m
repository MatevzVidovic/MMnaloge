

function X = presekPloskev(f1, gradf1, f2, gradf2, X0, h, n, tol = 1e-10, maxit = 100)
% f1 in f2 sta funkciji R^3 -> R.
% gradf1 in gradf2 sta funkciji pripadajocih gradientov. Vracata naj vrsticne vektorje
% (fino bi bilo, da tudi za stolpicne vektorje deluje)
% X0 je zacetni ugib, ki mora biti se popravljen, da postane x0
% h je velikost koraka (uporablja se v funkciji naslednjiUgib() )
% n je stevilo korakov
% tol in maxit sta omejitvi v Newtonovi metodi (uporablja se v funkciji popravekNewton() )
%
% funkcija presekPloskev poisce krivuljo, na kateri se ploskvi sekata.
% To stori tako, da se najprej za majhen korak premakne v smeri,
% v kateri nobena od funkcij ne spreminja vrednosti. (to ji poda funkcija smerniVektor() ) Tako dobi naslednji ugib.
% To opravi funkcija naslednjiUgib().
%
% Na tej tocki dobimo tudi vektor nadaljne smeri v ugibu. (s pomocjo funkcije smerniVektor())
%
% In potem pogledamo ravnino, ki je pravokotna na vektor smeri v ugibu, na kateri leži ugib.
% Na tej ravnini poiscemo tocko, kjer sta f1 in f2 enaki 0 (ta tocka lezi na obeh ploskvah)
% (Pogoj, da lezi na tej ploskvi, dobimo s h(x) = 0, kjer h(x) := v * x - v * y (y je ugib, v je vektor nadaljne smeri)
% Tako z Newtonovo metodo izboljšamo ta ugib.
% To stori funkcija popravekNewton()
%
% S ponavljanjem teh dveh korakov dobimo zaporedne tocke zeljene krivulje preseka ploskev.



  % definiranje G in JG
  % Ustvarimo najsirso konstrukcijo funkcije.
  % Potem zozamo na bolj uporabno konstrukcijo funkcije.
  % Ob vsaki uporabi je funkcija pred vnosom v newtonovo metodo ponovno konstruirana tako,
  % da je odvisna le od x (ostali parametri so podani)
  function z = konstrukcijaG(fone, ftwo, x, y, v)
    z(1, 1) = fone(x);
    z(2, 1) = ftwo(x);
    z(3, 1) = v' * x - v' * y;
  endfunction

  function z = konstrukcijaJG(gradfone, gradftwo, x, v)
    z(1, :) = gradfone(x);
    z(2, :) = gradftwo(x);
    z(3, :) = v';
  endfunction

  konstrG = @(x, y, v) (konstrukcijaG(f1, f2, x, y, v));
  konstrJG = @(x, v) (konstrukcijaJG(gradf1, gradf2, x, v));






  function v = konstrukcijaSmerniVektor(gradfone, gradftwo, x)
  % funkcija smerniVektor() iz funkcij gradientov izracuna vektor, v katerega smeri
  % se vrednost nobene od funkcij ne spreminja (pravokotno na oba gradienta v dani tocki)
    v0 = cross(gradfone(x)', gradftwo(x)');
    v = v0 ./ norm(v0);
  endfunction

  smerniVektor = @(x)(konstrukcijaSmerniVektor(gradf1, gradf2, x));




  % Definiranje tovrstne funkcije na pametnejsi nacin, kot sem prejsnje funkcije,
  % saj ni treba dodatno definirati funkcije s specificnimi parametri.
  % Slabsi primer je spodaj v komentarju.
  function y = naslednjiUgib(x, F = smerniVektor, step = h)
    % y je naslednji ugib po linearni aproksimaciji krivulje.
    y = x + step * F(x);
  endfunction
%
%  function y = konstrukcijaNaslednjiUgib(x, F, step)
%    y = x + step * F(x);
%  endfunction
%
%  naslednjiUgib = @(x) (konstrukcijaNaslednjiUgib(x, smerniVektor, h));






  function [X, n] = popravekNewton(F, JF, X0, tol = 1e-10, maxit = 100)
    %[X, n] = popravekNewton(F, JF, X0, tol, maxit) solves the (nonlinear)
    %system F(X) = 0 using the Newton's iteration with initial
    %guess X0. (JF is the Jacobi matrix of F.)

    for n = 1:maxit
      %Execute one step of Newton's iteration...
      X = X0 - feval(JF, X0)\feval(F, X0);
      %... and check if the new approximation is within prescribed tolerance.
      if(norm(X - X0) < tol)
        break;
      endif
      X0 = X;
    endfor

    %A warning in case the last approximation is not within specified tolerance.
    if(n == maxit)
      warning("no convergence after maxit iterations")
    endif

  endfunction















  % popravek prvotnega ugiba in zapis v izhodno matriko
  smerniVektor(X0);
  G = @(x) (konstrG(x, X0, smerniVektor(X0)));
  JG = @(x) (konstrJG(x, smerniVektor(X0)));

  X( : , 1) = popravekNewton(G, JG, X0, tol, maxit);


  % for zanka ugib, popravek, zapis
  for i = 1:n
    v0 = X( : , columns(X));
    v1 = naslednjiUgib(v0);
    G = @(x) (konstrG(x, v1, smerniVektor(v1)));
    JG = @(x) (konstrJG(x, smerniVektor(v1)));
    v2 = popravekNewton(G, JG, v1, tol, maxit);
    X( : , columns(X) + 1) = v2;
  endfor




%!test
%! F1 = @(X) (X(1) + X(2));
%! F2 = @(X) (-X(1) + X(2));
%! gradF1 = @(X) ([1, 1, 0]);
%! gradF2 = @(X) ([-1, 1, 0]);
%! Y = presekPloskev(F1, gradF1, F2, gradF2, [0.01; 0.02; 0.4], 0.1, 300, 1e-10, 100);
%! assert(Y(1:2, 1), zeros(2, 1) , 1e-10);
%! for i = 1:301
%!  assert(Y(1:2, i), zeros(2, 1) , 1e-10);
%! endfor
%!
%!

%!test
%! F1 = @(X) (X(1) + X(2) - 1);
%! F2 = @(X) (-X(1) + X(2) + 1);
%! gradF1 = @(X) ([1, 1, 0]);
%! gradF2 = @(X) ([-1, 1, 0]);
%! Y = presekPloskev(F1, gradF1, F2, gradF2, [0.01; 0.02; 0.4], 0.1, 300, 1e-10, 100);
%! assert(Y(1, 1), 1 , 1e-10);
%! assert(Y(2, 1), 0 , 1e-10);
%! for i = 1:301
%!  assert(Y(1, i), 1 , 1e-10);
%!  assert(Y(2, i), 0 , 1e-10);
%! endfor


%[0.8; 4.8; 2.2]
% n2 * r2 = 5
%!test
%! F1 = @(X) (X(1)^2 + X(3)^2 - X(2));
%! F2 = @(X) (X(2) - 5);
%! gradF1 = @(X) ([2*X(1), -1, 2*X(3)]);
%! gradF2 = @(X) ([0, 1, 0]);
%! Y = presekPloskev(F1, gradF1, F2, gradF2, [0.8; 0.8; 0.2], 0.1, 300, 1e-10, 100);
%! assert(Y(1, 1)^2 + Y(3, 1)^2, 5 , 1e-10);
%! assert(Y(2, 1), 5, 1e-10)
%! for i = 1:301
%!  assert(Y(1, i)^2 + Y(3, i)^2, 5 , 1e-10);
%!  assert(Y(2, i), 5, 1e-10)
%! endfor

%!test
%! F1 = @(X) ((X(1)/2)^2 + (X(3)/5)^2 - X(2));
%! F2 = @(X) (X(2) - 7);
%! gradF1 = @(X) ([X(1)/2, -1, 2*X(3)/25]);
%! gradF2 = @(X) ([0, 1, 0]);
%! Y = presekPloskev(F1, gradF1, F2, gradF2, [0.8; 0.8; 0.2], 0.1, 700, 1e-10, 100);
%! assert((Y(1, 1)/2)^2 + (Y(3, 1)/5)^2, 7 , 1e-10);
%! assert(Y(2, 1), 7, 1e-10)
%! for i = 1:701
%!  assert((Y(1, i)/2)^2 + (Y(3, i)/5)^2, 7 , 1e-10);
%!  assert(Y(2, i), 7, 1e-10)
%! endfor

%!test
%! F1 = @(X) ((X(1)/2)^2 + (X(3)/5)^2 - X(2));
%! F2 = @(X) (X(2) - 1);
%! gradF1 = @(X) ([X(1)/2, -1, 2*X(3)/25]);
%! gradF2 = @(X) ([0, 1, 0]);
%! Y = presekPloskev(F1, gradF1, F2, gradF2, [0.8; 0.8; 0.2], 0.1, 300, 1e-10, 100);
%! assert((Y(1, 1)/2)^2 + (Y(3, 1)/5)^2, 1 , 1e-10);
%! assert(Y(2, 1), 1, 1e-10)
%! for i = 1:301
%!  assert((Y(1, i)/2)^2 + (Y(3, i)/5)^2, 1 , 1e-10);
%!  assert(Y(2, i), 1, 1e-10)
%! endfor

%!test
%! F1 = @(X) ((X(1)/2)^2 + (X(3)/5)^2 - X(2));
%! F2 = @(X) (X(2) + X(3) - 5);
%! gradF1 = @(X) ([X(1)/2, -1, 2*X(3)/25]);
%! gradF2 = @(X) ([0, 1, 1]);
%! Y = presekPloskev(F1, gradF1, F2, gradF2, [0.8; 0.8; 0.2], 0.1, 1200, 1e-10, 100);
%! for i = 1:1201
%!  assert((Y(1, i)/2)^2 + (Y(3, i)/5)^2, -Y(3, i) + 5, 1e-10);
%! endfor







endfunction
