GHCOPTS=-O2 -threaded -rtsopts -XBangPatterns

all: symFinder symFinder-4d

symFinder: symFinder.hs HPositive.hs
	ghc ${GHCOPTS} symFinder.hs

symFinder-4d: symFinder-4d.hs HPositive.hs
	ghc ${GHCOPTS} symFinder-4d.hs

clean:
	rm symFinder symFinder-4d
	rm *.o *.hi
