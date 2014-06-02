import numpy
import lin
from lin import *
import random
import sys, time
import pp

# Technical, general functions
def chunks (l, n):
	return (l[i:i+n] for i in range (0,len(l),n))

def tuples (l, n):
	""" 
		Generates all tuples of length n from the list l
	"""
	if n == 1: return map(lambda x: [x], l)
	else:
		t = []
		for tup in tuples (l,n-1):
			for el in l:
				t.append([el]+tup)
		return t

def indicies_of_rank1_vectors (tensors): 
	"""
		Return the list of indices representing rank1 vectors
	"""
	return map (lambda x: [x], simple_tensors)

def indicies_of_rank2_vectors (tensors):
	"""
		Return the list of indices representing rank1 vectors
	"""
	t = []
	for a in tensors:
		for b in tensors:
			if a[0] != b[0] and a[1] != b[1]:
				t.append ([a,b])
	return t

def orthogonal_proj_subset (l, p, n):
	"""
		Returns the list of projectors of rank n that are orthogonal to p 
		and are build from elements of list l.
	"""
	if n == 0: return [p]
	else:
		ret = []
		for pp in l:
			if lin.is_orthogonal (pp,p): 
				ret.extend (orthogonal_proj_subset(l,pp+p,n-1))
		return ret


def generate_list_of_candidates (plist1,plist2,q,stop=-1,verbose=False):
	"""
		Function returns the list of potentially interesting
		e-symmetries. Projector e is formed from to projectors:
		one from list plist1 and the second from plist2. Resulting
		e must commute with projector q.

		Optional stop argument stops the function after checking 
		stop combinations.
	"""
	ret = []
	i = 0
	j = 0
	k = 0
	total = len(plist1)*len(plist2)
	for p1 in plist1:
		for p2 in plist2:
			i = i+1
			if verbose == True and i%10000 == 0: 
				print  "\r{0:6.2f}".format(i*100.0/total)+"% and found",\
					   k,"good projectors",j,"were promising",
				sys.stdout.flush()
			if i==stop: return ret
			if lin.is_orthogonal (p1,p2) and lin.does_commute (p1+p2,q):
				if lin.random_bp_test(p1+p2-2*q):
					ret.append (p1+p2-2*q)
					k = k+1
				j = j+1
	return ret

def generate_list_of_canditates_notbp (plist1,plist2,q,stop=-1):
	ret = []
	i = 0
	j = 0
	k = 0
	total = len(plist1)*len(plist2)
	for p1 in plist1:
		for p2 in plist2:
			i = i+1
			if i%100000 == 0: 
				print  i*100/total,"% and found",\
					   j,"good projectors"
			if i==stop: return ret
			if is_orthogonal (p1,p2) and does_commute (p1+p2,q):
				ret.append (p1+p2-2*q)
				j = j+1
	return ret

def find_good_esymmetries (n,q,stop=-1,erank=5):
	"""
		Returns the list of good e-symmetries such that e projector
		is made from n rank1 projectors and erank-n rank 2 projectors.
		Resulting e-symmetry is e - 2q.
	"""
	print "Generating projectors",
	sys.stdout.flush()

	start_time = time.time()

	opr1 = orthogonal_proj_subset (rank1_projectors, zero9, n)
	opr2 = orthogonal_proj_subset (rank2_projectors, zero9, erank-n)
	
	print "done in "+"{0:.2}".format(time.time() - start_time)+"s"

	print "Generating list of candidates "
	start_time = time.time()

	ret = generate_list_of_candidates (opr1,opr2,q,stop,verbose=True)
	print "done in "+"{0:.2f}".format(time.time() - start_time)+"s"

	return ret

def parallel_find_good_esymmetries (n,q,stop=-1,erank=5,chunksize=100,swapping=False):
	"""
		Returns the list of good e-symmetries such that e projector
		is made from n rank1 projectors and erank-n rank 2 projectors.
		Resulting e-symmetry is e - 2q.
	"""
	print "Generating projectors ",
	sys.stdout.flush()

	start_time = time.time()

	job_server = pp.Server (secret='none')
	
	start_time = time.time()
	
	job1 = job_server.submit (orthogonal_proj_subset,
							 (rank1_projectors, zero9, n),
							 (),
							 ("numpy","lin"))
	job2 = job_server.submit (orthogonal_proj_subset,
							 (rank2_projectors, zero9, erank-n),
							 (),
							 ("numpy","lin"))

	opr1 = job1()
	opr2 = job2()

	print "done in "+"{0:.2}".format(time.time() - start_time)+"s"
	print "Length of primary list:   ", len(opr1)
	print "Length of secondary list: ", len(opr2)

	if swapping:
		print "Swapping lists"
		tmp = opr1
		opr1 = opr2
		opr2 = tmp

	slices = chunks (opr2, chunksize)
	snum = len(opr2)/chunksize + (0 if len(opr2)%chunksize==0 else 1)
	
	total = len(opr1)*len(opr2)
	print "Okay, we've got"+" {0:,} ".format(total)+"combinations to check. Let's get to work!"
	print "Generating list of candidates "
	start_time = time.time()

	jobs = [job_server.submit(generate_list_of_candidates,(opr1,s,q,stop),(),("numpy","lin")) 
			for s in slices]

	ret = []
	i = 0
	k = 0
	for job in jobs:
		jr = job()
		k = k+len(jr)
		i = i+1
		print  "\r{0:6.2f}".format(i/(1.0*snum)*100.0)+"% and found",\
	  		    k,"good projectors",
	  	sys.stdout.flush()
		ret.extend (jr)

	print "\ndone in "+"{0:.2f}".format(time.time() - start_time)+"s"

	job_server.destroy()
	return ret

def filter_non_coCP (l):
	return [s for s in l if not lin.is_positive(partial_transpose(s))]

def filter_bp (l,n=10000):
	return [s for s in l if lin.random_bp_test (s,n)]

def parallel_filter_non_coCP (l, chunksize=100):
	slices = chunks (l, chunksize)
	snum = len(l)/chunksize + (0 if len(l)%chunksize==0 else 1)
	
	job_server = pp.Server (secret='none')
	jobs = [job_server.submit(filter_non_coCP,(s,),(partial_transpose,),("numpy","lin")) 
		for s in slices]
		
	ret = []
	i = 0
	for job in jobs:
		ret.extend (job())
		i = i+1
		print  "\r{0:6.2f}".format(i/(1.0*snum)*100.0)+"% done",
		sys.stdout.flush()

	job_server.destroy()
	return ret

def parallel_filter_bp (l,n=10000,chunksize=100):
	slices = chunks (l, chunksize)
	snum = len(l)/chunksize + (0 if len(l)%chunksize==0 else 1)
	
	job_server = pp.Server (secret='none')
	jobs = [job_server.submit(filter_bp,(s,n),(),("numpy","lin")) for s in slices]
		
	ret = []
	i = 0
	for job in jobs:
		ret.extend (job())
		i = i+1
		print  "\r{0:6.2f}".format(i/(1.0*snum)*100.0)+"% done",
		sys.stdout.flush()
		
	job_server.destroy()
	return ret
		
# Random number generator
random.seed()

#
# Few important matrices
#

# Choi matrix of transposition
wsym = numpy.array([
				[1,0,0,0,0,0,0,0,0], 
				[0,0,0,1,0,0,0,0,0],
				[0,0,0,0,0,0,1,0,0],
				[0,1,0,0,0,0,0,0,0],
				[0,0,0,0,1,0,0,0,0],
				[0,0,0,0,0,0,0,1,0],
				[0,0,1,0,0,0,0,0,0],
				[0,0,0,0,0,1,0,0,0],
				[0,0,0,0,0,0,0,0,1]], dtype=float)

# Projector on maximally entangled vector Schmidt ra nk 3
q3 = 1.0/3.0*projector(sum(map(make_vector_from_indices, [[0,0],[1,1],[2,2]])))
# Projector on entangled vector Schmidt rank 2
q2 = 0.5*projector(sum(map(make_vector_from_indices, [[0,0],[1,1]])))
## 2D projector
qq2 = 0.5*projector(sum(map(make_vector_from_indices, [[0,0],[1,1]]))) + \
		0.5*projector(make_vector_from_indices([0,2])-make_vector_from_indices([2,0]))

qq2a = 0.5*projector(make_vector_from_indices([0,0])+make_vector_from_indices([1,1])) + \
		0.5*projector(make_vector_from_indices([0,1])-make_vector_from_indices([1,0]))
		
qq2b = 0.5*projector(make_vector_from_indices([0,0])+make_vector_from_indices([1,1])) + \
		0.5*projector(make_vector_from_indices([1,2])+make_vector_from_indices([2,1]))
		
qq2c = 0.5*projector(make_vector_from_indices([0,0])+make_vector_from_indices([1,1])) + \
		0.5*projector(make_vector_from_indices([1,2])-make_vector_from_indices([2,0]))

qq2d = 0.5*projector(make_vector_from_indices([0,0])+make_vector_from_indices([1,1])) + \
		0.5*projector(make_vector_from_indices([0,2])+make_vector_from_indices([2,1]))

		
# Tuples of simple tensors
simple_tensors = tuples (range(3),2)
# List of indicies of all rank1, rank2 interesting vectors
rank1_vectors = indicies_of_rank1_vectors (simple_tensors)
rank2_vectors = indicies_of_rank2_vectors (simple_tensors)

# Projectors made from rank1 vectors lists
rank1_projectors = make_rank1_projectors (rank1_vectors)
rank2_projectors = make_rank2_projectors (rank2_vectors)

if __name__ == "__main__":
        # add the directory with ppworker.py to the path
        sys.path.append(os.path.dirname(__file__))
        sys.path.pop(0)
        f.close()
        wp = _WorkerProcess()
        wp.run()
        
