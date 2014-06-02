Supplement for PhD thesis
=============

Content of this repository is a supplementary
material for my PhD thesis "Structure of the set of positive
maps between low dimensional matrix algebras".
These are mainly source codes used during the research.

## Description

The file `thesis notebok.nb` is a Mathematica notebook
that contains symbolic calculations to complex to be
put explicitly in the thesis. It also illustrates how to
use `positive.m` library, which was created during the
research to simplify tasks. The source code of this
library is in `positive.nb` notebook, along with some
basic documentation. 

I don't attach original notebooks, as they often miss
the comments and thus would be unreadable. 

`HPositive.hs` is a Haskell module for brute-force search
of symmetries of the form Eq. (2.12) from the thesis.
Actual 'workers' are in `symFinder-4d.hs` for n = dim H = 4,
and `symFinder.hs` for experimenting with other settings.
Complexity of bruteforce search scales so fast, that
for n > 4 it is practically useless. I started developing
optimized version which could give some results for bigger
n, but it is still not completed.

Results of the `symFinder-4d.hs` for n = 4 are uploaded
in `candidates-*.m` files, in the format that can be directly
imported to the Mathematica.

Finally, python sources have been used to find any
interesting example of partial symmetry in n = 3 case.
It did not succeeded. The code was written quite a long
time ago and is rather poorly documented, but I submit it,
as I reference results obtained by it.

## Installation

Mathematica notebooks doesn't need any installations.
Should be simply opened by the recent version of Mathematica.
The only thing is that you shouldn't separate `positive.m` from
the `thesis notebook.nb` as the latter loads the former, assuming
that it is in the same directory.

For Haskell code, you need to install [Haskell](www.haskell.org) platform.
We refer to that webpage for detailed install instructions,
depending on the platform you use.

You will also need to use custom version of `hmatrix` Haskell library.
This can be cloned from git repository:

```sh
$ git clone https://github.com/ttylec/hmatrix.git
```

then from within `hmatrix` directory run

```sh
$ cabal configure
$ cabal build
$ cabal install
```

and everything should work.

To compile the source code from this repository simply run `make`.

Python source code does not require any installation steps
(you might need to install parallel python module).
