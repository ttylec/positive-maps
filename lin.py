import numpy
import random

def hermitian_conjugate (a): return numpy.transpose(numpy.conjugate(a))

def basis_vector (i,dim=3):
	e = numpy.zeros ( (dim,1) )
	e[i,0] = 1.0
	return e

def projector (x):
	"""
		Return projector on vector x
	"""
	return numpy.dot ( x, hermitian_conjugate(x) )

def inner_product (x,y): return numpy.inner (numpy.conjugate(x),y)

def is_orthogonal (p,q): return numpy.array_equal(numpy.dot(p,q), zero9)
	
def does_commute (a,b): return numpy.array_equal(numpy.dot(a,b)-numpy.dot(b,a),zero9)

def is_positive (a): return min(numpy.linalg.eigvalsh(a))>-0.001

def random_vector (dim=3):
	"""
		Generate random vector of dimension dim
	"""
	return numpy.array(map (lambda x: random.random()+random.random()*1.0j, range(dim)))

def random_bp_test (a, n=1000):
	"""
		Random block-positivity test of matrix a, with n trials
	"""
	i = n
	while (n>0):
		x = random_vector ()
		y = random_vector ()
		xy = numpy.kron(x,y)
		bpcond = inner_product (xy, numpy.dot(a,xy))
		if bpcond.real<-0.001: return False
		n = n-1
	return True

def make_vector_from_indices (i):
	"""
		Returns vector e_i\otimes e_j for pair of indicies (i,j)
	"""
	return numpy.kron(basis_vector(i[0]),basis_vector(i[1]))

def make_rank1_projectors (r2v):
	"""
		Generate list of rank 1 projectors from list of indicies
	"""
	ret = []
	for v in r2v:
		x = make_vector_from_indices(v[0])
		ret.append(projector(x))
	return ret

def make_rank2_projectors (r2v):
	"""
		Generate list of rank 1 projectors from list of indicies	
	"""
	ret = []
	for v in r2v:
		x = make_vector_from_indices(v[0])
		y = make_vector_from_indices(v[1])
		ret.append(0.5*projector(x+y))
		ret.append(0.5*projector(x-y))
	return ret

def partial_transpose (a):
	ret = numpy.zeros((9,9));
	for i in range(3):
		for j in range(3):
			ret[3*i:3*(i+1),3*j:3*(j+1)] = numpy.transpose(a[3*i:3*(i+1),3*j:3*(j+1)])
	return ret

# zero 9x9 matrix
zero9 = numpy.zeros((9,9))
