#   Copyright (C) 2010   Marcus P da Silva http://github.com/marcusps
# 
#   License: Distributed under GPL 2.0
#            http://creativecommons.org/licenses/GPL/2.0/
#            http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.
import numpy as np

def ket( s ) :
  """
  This function returns a column vector numpy.matrix with entried given by nlist.
  """
  return bra( s ).T;

def bra( s ) :
  """
  This function returns a row vector numpy.matrix with entries given by nlist.
  """
  if s.__class__ == str:
    res = np.zeros( 2**len(s) );
    res[ int( s, 2 ) ] = 1;
    return np.matrix( res );
  else:
    return np.matrix( s );

def commutator( a, b ):
  """
  This function takes two matrices (of type numpy.matrix) and returns
  their commutator [a,b], defined as a*b-b*a.
  """
  return a*b - b*a;

