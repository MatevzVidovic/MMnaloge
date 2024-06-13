## Copyright (C) 2022 Uporabnik
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} zanimivosti (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Uporabnik <Uporabnik@MATEVZ-W10>
## Created: 2022-05-12

function retval = zanimivosti (input1, input2)





  % za zrihtat


   F = @(X) (X.^3 - sin(X));
   [P, Q] = presekKrivulj(@elipsaEna, @dotElipsaEna, [0, 2*pi], @sinus, @dotSinus, [0, 2*pi], 0.05);
   assert(F(P(1, : )), zeros(1, columns(P)) , 1e-10)



  function R = elipsaEna (inputVect)
    R(1, : ) = 3 * cos(inputVect);
    R(2, : ) = 1.5 * sin(inputVect);
  endfunction

  function R = dotElipsaEna (inputVect)
    R(1, : ) = -3 * sin(inputVect);
    R(2, : ) = 1.5 * cos(inputVect);
  endfunction

endfunction
