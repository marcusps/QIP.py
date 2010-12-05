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

from numpy import reshape
from scipy.sparse import kron, eye

def vec(m):
    ''' 
    This function returns the column-major Liouville
    representation of a linear operator.
    '''
    return reshape(m.T,(m.shape[0]*m.shape[1],1));
#end def  

#def vech(m):
#    return vec(m);
##end def

def vecinv( m, shape ):
    ''' 
    This function returns the natural
    representation of a linear operator with a certain shape, given
    its Liouville representation.
    '''
    return reshape(m,shape).T;
#end def

#def vechinv( m, shape ):
#    return vecinv(m,shape);
##end def

def liou(l, r):
    ''' 
    This function returns the column-major Liouville
    representation of a linear operator acting as L on the left and as
    R on the right.
    '''
    return kron(r.T, l)
#end def

def kraus2liou( l ):
    ''' 
    This function returns the column-major Liouville
    representation of a set of Kraus operators.
    '''
    if l.__class__ == [].__class__:
        return reduce(lambda x,y: x+y,map(lambda x: liou(x,x.H),l))
    else:
        return liou( l, l.H );
    #end if
#end def

def choi2liou( c ):
    ''' 
    This function returns the Choi representation of
    a column-major Liouville representation of a superoperator.
    '''
    _s = c.shape;
    _M = c;

    _index = lambda a,b: (a-1)*sqrt(s[0])+b;
    
    for n in range( sqrt(s[0]) ):
        for m in range( sqrt(s[0]) ):
            for o in range( sqrt(s[0]) ):
                for p in range( n+1, sqrt(s[0]) ):
                    # this swaps M_{nm,op} with M_{pm,on}
                    temp = M[ index(n,m), index(o,p) ];
                    M[ index(n,m), index(o,p) ] = M[ index(p,m), index(o,n) ];
                    M[ index(p,m), index(o,n) ] = temp;
                    
##end def


def dissipator(a):
    ''' 
    This function returns the column-major Liouville
    representation of a dissipator superoperator.
    '''
    idm = eye(a.shape[0])
    return liou(a,a.H)-0.5*liou(a.H*a,idm)-0.5*liou(idm,a.H*a)
#end def

def hamiltonian(h):
    ''' 
    This function returns the column-major Liouville
    representation of a Hamiltonian superoperator ( the commutator of an operator
    with the Hamiltonian in the Schroedinger picture ).
    '''
    idm = eye(h.shape[0])
    return liou(h,idm)-liou(idm,h)
#end def

def dual(m):
    ''' 
    This function returns the dual of a superoperator in the Liouville representation.
    '''
    return m.H;
#end def
