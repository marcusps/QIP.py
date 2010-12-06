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

id = np.matrix([[1,0],[0,1]]);
x  = np.matrix([[0,1],[1,0]]);
y  = np.matrix([[0,-1j],[1j,0]]);
z  = np.matrix([[1,0],[0,-1]]);

def _char2pauli( c ):
  """ 
  Convert a character to the corresponding Pauli operator. If the
  character is not 'x', 'y', 'z' or 'i', an exception is raised.
  """
  if c == 'x' or c == 'X':
    return x
  elif c == 'y' or c == 'Y':
    return y
  elif c == 'z' or c == 'Z': 
    return z
  elif c == 'i' or c == 'I': 
    return id
  else:
    raise ValueError, ( 'Unrecognized Pauli operator \'' + c + '\'')

def _num2pauli( c ):
  """ 
  Convert a number to the corresponding Pauli operator. If the
  number is not 0, 1, 2 or 3, an exception is raised.
  """
  if c == 1:
    return x
  elif c == 3:
    return y
  elif c == 2: 
    return z
  elif c == 0: 
    return id
  else:
    raise ValueError, ( 'Unrecognized Pauli operator \'' + c + '\'' )

def string2pauli( s ):
  """
  Convert a string made up of the characters 'x', 'y', 'z' and 'i' (or
  their uppercase forms) into the matrix resulting from the tensor
  prodict (Kronecker product) of the corresponding Pauli matrices.
  """
  return reduce(np.kron,map(_char2pauli,s))

def list2pauli( s ):
  """
  Convert a list made up of the integers 0, 1, 2 and 3 into the matrix
  resulting from the tensor prodict (Kronecker product) of the
  corresponding Pauli matrices.
  """
  return reduce(np.kron,map(_num2pauli,s))

