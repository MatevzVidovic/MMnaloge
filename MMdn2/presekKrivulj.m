

function [P, Q] = presekKrivulj (p, pdot, intp, q, qdot, intq, h)
  %funkcija presekKrivulj vzame vektorski funkciji p in q, ki slikata iz R v R^2
  %ter vektorski funkciji pdot in qdot, ki sta njuna odvoda.
  %Vzame tudi dvomestna vektorja intp ter intq, katerih prva komponenta je začetek intervala,
  %na katerem funkcija tece, druga komponenta pa konec intervala.
  % Prejme tudi stevilo h, ki predstavlja razmike v zacetni linearni aproksimaciji krivulj.

  %funkcija presekKrivulj vraca preseke prejetih krivulj v R^2
  %to stori tako, da najprej za vsako od krivulj ustvari linearno aproksimacijo te funkcije.
  %To linearno aproksimacijo dobi tako, da interval funkcije razreze na dele sirine h.
  %Potem za vsako od komponentnih funkcij izracuna vrednost za vsako od delilnih tock.
  %In potem za vsakega od dobljenih odsekov funkcije p preveri, ce se seka s katerimkoli od odsekov funkcije q.
  %Kjer se, pridobimo koordinate presecisca teh linearnih aproksimacij.
  %Potem za vsako pridobljeno presecisce izvede Newtonovo metodo, kjer prej pridobljene koordinate uporabi za prvotno ugibanje.





  %ustvarjanje delitve na intervalu

  numP = ceil((intp(2) - intp(1)) / h)+1;
  %ceil, ker bolje vec, kot pa premalo. +1 pa, ker moramo upostevati zacetek oziroma konec kot dodatno delilno tocko.
  intVectP = linspace(intp(1), intp(2), numP);

  %delni test
  % intp = [1; 10.5];
  % h = 1;

  numQ = ceil((intq(2) - intq(1)) / h)+1;
  intVectQ = linspace(intq(1), intq(2), numQ);



  %to bosta matriki z dvema vrsticama - prva za x in druga za y
  linP = p(intVectP);
  linQ = q(intVectQ);




  % ima 6 vrstic. Prva vrstica je x koordinata presecisca, druga je y.
  % tretja je i v tej iteraciji (aproksimacija t za funkcijo p)
  % cetrta je j v tej iteraciji (aproksimacija t za funkcijo q)
  % peta je parameter presecisca za p
  % sesta je parameter presecisca za q
  linPres = [];



  %preverjanje sekanja
  for i = 1:(length(linP)-1)
    for j = 1:(length(linQ)-1)

      v(1, 1) = linQ(1, j) - linP(1, i);
      v(2, 1) = linQ(2, j) - linP(2, i);

      M(1,1) = linP(1, i+1) - linP(1, i);
      M(2,1) = linP(2, i+1) - linP(2, i);

      M(1,2) = -(linQ(1, j+1) - linQ(1, j));
      M(2,2) = -(linQ(2, j+1) - linQ(2, j));

      %test na tej tocki:
      % presekKrivulj(@kub, @dotKub, [-1, 5], @kvadrat, @dotKvadrat, [-1, 5], 1);
      % pa samo odkomentiraj naslednji line:
      %M


      if (M(:,1) == M(:,2))

        % tu bi moral implementirat special cases, ampak jih bom pac izpustil
        % saj numericno gledano prakticno vedno, ko sta vekotrja enaka, ne bo sekanja

      else

        koef = M\v;

        if (0 <= koef(1) && koef(1) < 1 && 0 <= koef(2) && koef(2) < 1)
          linPres(1, columns(linPres)+1) = linP(1, i) + koef(1) * (linP(1, i+1) - linP(1, i));
          linPres(2, columns(linPres)) = linP(2, i) + koef(1) * (linP(2, i+1) - linP(2, i));
          linPres(3, columns(linPres)) = i;
          linPres(4, columns(linPres)) = j;
          linPres(5, columns(linPres)) = (intp(1) + (i-1) * h + koef(1) * h);
          linPres(6, columns(linPres)) = (intq(1) + (j-1) * h + koef(2) * h);

        endif

      endif

    endfor
  endfor

  Q = linPres(1:2, : );

  %izris grafa
  figure('Name','Linearna aproksimacija');
  plot(linP(1, :), linP(2, : ))
  hold on
  plot(linQ(1, : ), linQ(2, : ))
  plot(linPres(1, :), linPres(2, : ), "or")
  axis equal






  %delni testi
  %tukaj je problem, ker je odvod pri obeh 0 v tocki 0, pa mogoce se kaj
  % presekKrivulj(@kub, @dotKub, [-1, 2], @kvadrat, @dotKvadrat, [-1, 2], 0.1);
  % presekKrivulj(@kub, @dotKub, [-1, 2], @kvadrat, @dotKvadrat, [-1, 2], 0.01);

  % presekKrivulj(@kub, @dotKub, [-1, 2], @sinus, @dotSinus, [-1, 2], 0.01);
  % presekKrivulj(@kub, @dotKub, [0.5, 1.5], @sinus, @dotSinus, [0.5, 1.5], 0.01);
  % presekKrivulj(@cosinus, @dotCosinus, [-2*pi,2*pi], @sinus, @dotSinus, [-2*pi, 2*pi], 0.3);
  % presekKrivulj(@kub, @dotKub, [-1, 2], @sinus, @dotSinus, [-1, 2], 0.1);


















  function [X, n, correct] = newton(F, JF, X0, tol = 1e-10, maxit = 100)
  %[X, n, correct] = newton(F, JF, X0, tol, maxit) solves the (nonlinear)
  %system F(X) = 0 using the Newton's iteration with initial
  %guess X0. (JF is the Jacobi matrix of F.)
  % n is the number of steps, and correct is boolean for if the method even converged

  correct = 0;

  for n = 1:maxit
	  %Execute one step of Newton's iteration...
	  X = X0 - feval(JF, X0)\feval(F, X0);
	  %... and check if the new approximation is within prescribed tolerance.
	  if(norm(X - X0) < tol)
      correct = 1;
		  break;
	  endif
	  X0 = X;
   endfor

    %A warning in case the last approximation is not within specified tolerance.
    if(n == maxit)
      warning("no convergence after maxit iterations")
    endif

  endfunction





  % izboljsava z newtonovo metodo

  F = @(X) (p(X(1))-q(X(2)));

  function [M] = jacobi(X, fone, ftwo)
    M(:, 1) = fone(X(1));
    M(:, 2) = -(ftwo(X(2)));
  endfunction

  JF = @(X) jacobi(X, pdot, qdot);


  izboljsanLinPres = [];
  for i = 1:columns(linPres)
    zacX = linPres(5:6, i);
    [V, stPon, isCorrect] = newton(F, JF, zacX);
    if (isCorrect)
      izboljsanLinPres(:, columns(izboljsanLinPres)+1) = V;
    endif
  endfor




  % pretvorba parametrov presecisc v koordinate
  P = [];
  for i = 1:columns(izboljsanLinPres)
    P(:, columns(P)+1) = p(izboljsanLinPres(1, i));
  endfor



  %izris grafa
  % fplot je bil tezaven za izris krivulj, saj je bolj namenjen izrisu funkcij R -> R in R^2 - R.
  % Stalno je zelel prikazati tudi interval, od koder se slika, kar je izpadlo kot dodaten izris simetrale lihih kvadrantov.
  % Iz tega razloga sem uporabil kar plot ze prej napisane linearne aproksimacije, saj naredi prakticno isto kot fplot,
  % le da to naredi bolje, saj ga ne mede dejstvo, da je funkcija vektorska.

  figure('Name','Newtonova metoda');
  %fplot(p, intp);
  plot(linP(1, :), linP(2, : ))
  hold on
  %fplot(q, intq);
  plot(linQ(1, : ), linQ(2, : ))
  plot(P(1, :), P(2, : ), "or")
  axis equal








  %testiranje:

  function R = kub (inputVect)
    R(1, : ) = inputVect;
    R(2, : ) = inputVect.^3;
  endfunction

  function R = dotKub (inputVect)
    R(1, : ) = inputVect;
    R(2, : ) = (inputVect.^2).*3;
  endfunction

  function R = sinus (inputVect)
    R(1, : ) = inputVect;
    R(2, : ) = sin(inputVect);
  endfunction

  function R = dotSinus (inputVect)
    R(1, : ) = inputVect;
    R(2, : ) = cos(inputVect);
  endfunction

  function R = cosinus (inputVect)
    R(1, : ) = inputVect;
    R(2, : ) = cos(inputVect);
  endfunction

  function R = dotCosinus (inputVect)
    R(1, : ) = inputVect;
    R(2, : ) = -sin(inputVect);
  endfunction





%!test
%! F = @(X) (cos(X) - sin(X));
%! [P, Q] = presekKrivulj(@cosinus, @dotCosinus, [-2*pi, 2*pi], @sinus, @dotSinus, [-2*pi, 2*pi], 0.1);
%! assert(F(P(1, : )), zeros(1, columns(P)) , 1e-10);

%!test
%! F = @(X) (X.^3 - cos(X));
%! [P, Q] = presekKrivulj(@cosinus, @dotCosinus, [-1, 1.5], @kub, @dotKub, [-1, 1.5], 0.1);
%! assert(F(P(1, : )), zeros(1, columns(P)) , 1e-10);

%!test
%! F = @(X) (X.^3 - sin(X));
%! [P, Q] = presekKrivulj(@kub, @dotKub, [-1, 2], @sinus, @dotSinus, [-1, 2], 0.02);
%! assert(F(P(1, : )), zeros(1, columns(P)) , 1e-10);




endfunction
