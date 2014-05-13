module HPositive (
                 basis, tensorBasis, tensor, proj,
                 basicSchmidtRankSet,
                 schmidtRank, randomBpTest,
                 randomBpTest3, randomBpTest4, randomBpTest5, randomBpTest6, randomBpTest7,
                 symmetryWorker,
                 findSymmetries,
                 parallelConstructSymmetry,
                 parallelConstructSymmetryStepII,
                 constructSymmetry,
                 -- findSymmetries',
                 addCounts,
                 phases2R,
                 phases2C,
                 phases3R
                 ) where

import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.LAPACK
import System.Random
import Control.Parallel.Strategies
import Control.DeepSeq
import Data.List
import Text.Printf
import System.IO

import qualified Debug.Trace as Debug
debug = flip Debug.trace

-- old version of hmatrix
cdot = (<.>)

--
-- Basic linear algebra stuff
-- 
type HVec = Vector (Complex Double)
type Traces = Vector (Complex Double)
type Operator = Matrix (Complex Double)
type TList = [Complex Double]

-- list of basis vectors in tensor product of two dim-dimensional Hilbert
-- spaces
tensorBasis :: Int -> [HVec]
tensorBasis dim = [(e i) `tensor` (e j) | i <- [1..dim], j <- [1..dim]]
    where
        e i = basis dim i

-- tensor product of two vectors
tensor :: HVec -> HVec -> HVec
tensor x y = flatten $ kronecker (asColumn x) (asColumn y)

-- i-th basis vector in dim-dimensional Hilbert space 
basis :: Int -> Int-> HVec
basis dim i = fromList $ prefix ++ [unit] ++ suffix
    where
        unit = 1.0 :: Complex Double
        prefix = replicate (i-1) 0.0 :: [Complex Double]
        suffix = replicate (dim-i) 0.0 :: [Complex Double]

-- Schmidt rank of vector
schmidtRank :: HVec -> Int
schmidtRank v = matrixRank coeff_mat
    where
        d = (floor . sqrt . fromIntegral . dim) v
        coeff_mat = buildMatrix d d coeff
        coeff (i, j) = cdot eij v
            where
                eij = tensor (basis d (i+1)) (basis d (j+1))

matrixRank :: Operator -> Int
matrixRank = length . filter (>accuracy) . toList . svC

trace :: Operator -> Complex Double
trace a = sum $ toList $ takeDiag a 

-- projector onto given vector (normalized)
proj :: HVec -> Operator
proj v = scale nn $ cv <> (ctrans cv)
    where 
        cv = asColumn v
        nn = 1/(norm2 v)^2 :+ 0.0

-- return True if given list of vectors is linearly independent
isIndependentSubset :: [HVec] -> Bool
isIndependentSubset vs = length vs == (matrixRank . fromRows) vs

-- return True if given list of vectors is mutually orthogonal
isOrthoSubset :: [HVec] -> Bool
isOrthoSubset []   = True
isOrthoSubset (v:[]) = True
isOrthoSubset (v:vs) = orthoToAll && isOrthoSubset vs 
    where
        orthoToAll = all equalsZero $ map (cdot v) vs

-- random test of block positivity
randomBpTest :: Int -> Operator -> Bool
randomBpTest trials a = all (>= 0) inner
    where
        inner = map (\v -> realPart (v `cdot` (a <> v))) vlist
        d = (floor . sqrt . fromIntegral . rows) a
        vlist = map (\s -> randomTensorVector s d) [1..trials]

randomBpTest3 :: Int -> Operator -> Bool
randomBpTest3 trials a = all (>= 0) inner
    where
        inner = map (\v -> realPart (v `cdot` (a <> v))) $ take trials randomVectorPool3

randomBpTest4 :: Int -> Operator -> Bool
randomBpTest4 trials a = all (>= 0) inner
    where
        inner = map (\v -> realPart (v `cdot` (a <> v))) $ take trials randomVectorPool4

randomBpTest5 :: Int -> Operator -> Bool
randomBpTest5 trials a = all (>= 0) inner
    where
        inner = map (\v -> realPart (v `cdot` (a <> v))) $ take trials randomVectorPool5

randomBpTest6 :: Int -> Operator -> Bool
randomBpTest6 trials a = all (>= 0) inner
    where
        inner = map (\v -> realPart (v `cdot` (a <> v))) $ take trials randomVectorPool6

randomBpTest7 :: Int -> Operator -> Bool
randomBpTest7 trials a = all (>= 0) inner
    where
        inner = map (\v -> realPart (v `cdot` (a <> v))) $ take trials randomVectorPool7

randomComplexVector :: Int -> Int -> HVec
randomComplexVector seed d = zipVectorWith (:+) realPart imaginaryPart
    where
        [realPart, imaginaryPart] = takesV [d, d] $ mapVector (\x -> (1-2*x)) $ randomVector seed Uniform (2*d)

randomTensorVector :: Int -> Int -> HVec
randomTensorVector seed d = tensor x y
    where
        [x, y] = takesV [d, d] $ randomComplexVector seed (2*d)

-- We optimize bp checks by defining global pool of random vectors
randomVectorPool3 = map (\s -> randomTensorVector s 3) [1..]
randomVectorPool4 = map (\s -> randomTensorVector s 4) [1..]
randomVectorPool5 = map (\s -> randomTensorVector s 5) [1..]
randomVectorPool6 = map (\s -> randomTensorVector s 6) [1..]
randomVectorPool7 = map (\s -> randomTensorVector s 7) [1..]

--
-- Constructing symmetries
--

--
-- Unital method
--
-- This should be very fast but will work only for fixed ranks
--

-- This is outer loop for bruteforce check
-- Outside loop is called sequentially, to give some feedback about
-- the progress
parallelConstructSymmetry :: (Operator -> Bool) -> Int -> [HVec] -> IO ()
parallelConstructSymmetry bptester aim (x:xs) = do
        let
            q0 = proj x
            counts@(c0:cs) = snd . head . addCounts $ [x]
            qq = (q0, counts)
            xcs = addCounts xs
            oxs = xcs `onlyOrthoToProj'` q0
            goodxs = filter (\(x, c) -> (head c) > 0) $ oxs
            nextxs = map (\(x, c) -> (x, tail c)) $ filter (\(x, c) -> (head c) == 0) $ oxs
            qs = map (expandByProj qq) $ filter (isGoodMix aim qq) ps
            ps = map makeQ $ filter (isOrthoSubset . fst. unzip) $ subsets (aim - c0) goodxs
        hPutStrLn stderr $ printf "%d to check." $ length qs
        mapM_ (\(q, idx) -> parallelConstructSymmetryStepII idx bptester aim q nextxs) $ zip qs [1..]

-- This schedules search of symmetries and stores the result in
-- file if anything was found.
-- This is a second step of search: for given 'preconstructed' q it
-- executes construction of possible q's in parallel.
parallelConstructSymmetryStepII idx bptester aim qq@(q0, counts@(c0:cs)) xs = do
        let
            qs = map (expandByProj qq) $ filter (isGoodMix aim qq) ps
            ps = map makeQ $ filter (isOrthoSubset . fst. unzip) $ subsets (aim - c0) goodxs
            oxs = xs `onlyOrthoToProj'` q0
            goodxs = filter (\(x, c) -> (head c) > 0) $ oxs
            nextxs = map (\(x, c) -> (x, tail c)) $ filter (\(x, c) -> (head c) == 0) $ oxs
            sc = map (\q -> constructSymmetry aim q nextxs) qs `using` (parList rdeepseq)
            ss = filter bptester $ concat $ sc 
            fn = printf "candidates-%d.m" (idx :: Int) 
        hPutStrLn stderr $ printf "%d: need to check %d extensions" (idx :: Int) (length qs)
        if length ss > 0 
        then 
            do
                outf <- openFile fn WriteMode
                hPutStrLn outf $ printf "{"
                hPutStrLn outf $ intercalate ",\n" $ map mfMatrix $ ss
                hPutStrLn outf "}"
                hClose outf
        else
            do
                return ()
        hPutStrLn stderr $ printf "%d: found %d" (idx :: Int) (length ss)


expandByProj (q0, counts) (p, cs) = (q0 + p, tail $ sumLists counts cs)

constructSymmetry :: Int -> (Operator, [Int]) -> [(HVec, [Int])] -> [Operator]
constructSymmetry _ (q, []) _ = [makeS q] -- `debug` "all aims done: constructing symmetry"
constructSymmetry aim (q, cs) []
    | all (== aim) cs = [makeS q] -- `debug` "empty vector list, but aims fulfilled: constructing symmetry" 
    | otherwise = [] -- `debug` "empty vector list, but aims not fulfilled"
constructSymmetry aim qq@(q0, counts@(c0:cs)) xs
        | c0 == aim = constructSymmetry aim (q0, cs) nextxs -- `debug` ("aim fulfilled, skipping (" ++ ((show . length) cs) ++ " to go)")
        | qs == [] && all (== aim) counts = [makeS q0] -- `debug` "cannot extend, but aims fulfilled: constructing symmetry"
        | qs == [] && any (/= aim) counts = [] -- `debug` "cannot extend and aim is not fulfilled"
        | otherwise =  concat $ map (\q -> constructSymmetry aim q nextxs) $ qs --`debug` formatDbgMsg
        where
            qs = map (expandByProj qq) $ filter (isGoodMix aim qq) ps
            ps = map makeQ $ filter (isOrthoSubset . fst. unzip) $ subsets (aim - c0) goodxs
            oxs = xs `onlyOrthoToProj'` q0
            goodxs = filter (\(x, c) -> (head c) > 0) $ oxs
            nextxs = map (\(x, c) -> (x, tail c)) $ filter (\(x, c) -> (head c) == 0) $ oxs
            traceConstruction :: a -> a
            traceConstruction  
                | length counts == dim = Debug.trace (((show . length) qs) ++ " to check at level 1...")
                | length counts == dim - 1 = Debug.trace (((show . length) qs) ++ " to check at level 2...")
                | length counts == dim - 2 && dim > 4 = Debug.trace "recursive call at level 3..."
                | otherwise = id
            dim = round . sqrt . fromIntegral . rows $ q0

isGoodMix aim (_, counts) (_, cs) = (maximum sums <= aim) && (head sums == aim)
    where
        sums = sumLists cs counts

makeQ :: [(HVec, [Int])] -> (Operator, [Int])
makeQ vts = (foldl1 (+) $ map proj vs, foldl1 (sumLists) counts)
    where
        (vs, counts) = unzip vts

addCounts :: [HVec] -> [(HVec, [Int])]
addCounts vs = map (\x -> (x,  counts x)) $ vs
    where
        n = round . sqrt . fromIntegral . dim $ head vs
        counts x = map ((partialTraceCheck x) . (basis n)) $ [1..n]
        partialTraceCheck x y
            | (magnitude $ unitalTrace (proj x) (proj y)) > 0 = 1
            | otherwise = 0

-- Interface
findSymmetries bptester rank (x:xs) = do
        let
            q0 = proj x
            qs = map ((q0 + ) . proj) $ xs `onlyOrthoToProj` q0
        hPutStrLn stderr $ printf "%d will be checked" (length qs)
        mapM_ (symmetryWorker bptester xs rank) $ zip qs [1..]

symmetryWorker bptester xs rank (q, idx) = do
        let ssl = map (symmetryCandidates bptester xs (rank-3)) qs `using` (parList rdeepseq)
            qs = map ((q + ) . proj) $ xs `onlyOrthoToProj` q
            ss = ssl `seq` concat ssl
            fn = printf "candidates-%d.m" (idx :: Int) 
        outf <- openFile fn WriteMode
        hPutStrLn outf $ printf "{"
        hPutStrLn outf $ intercalate ",\n" $ map mfMatrix $ ss
        hPutStrLn outf "}"
        hPutStrLn stderr $ printf "%d checking, found %d" (idx :: Int) (length ss)

-- findSymmetries' bptester rank (x:xs) = do
--         let
--             xst = map (\x -> (x, unitalTraces x)) $ xs
--             q0 = (proj x, unitalTraces x)
--             qs = map (extendByVector q0) $ onlyPromisingVectors xst q0
--         hPutStrLn stderr $ printf "%d will be checked" (length qs)
--         mapM_ (symmetryWorker' bptester xst rank) $ zip qs [1..]
-- 
-- symmetryWorker' bptester xs rank (q, idx) = do
--         let ssl = map (symmetryCandidates' bptester xs (rank-3)) qs `using` (parList rdeepseq)
--             qs = map (extendByVector q) $ onlyPromisingVectors xs q
--             ss = ssl `seq` concat ssl
--             fn = printf "candidates-%d.m" (idx :: Int) 
--         outf <- openFile fn WriteMode
--         hPutStrLn outf $ printf "{"
--         hPutStrLn outf $ intercalate ",\n" $ map mfMatrix $ ss
--         hPutStrLn outf "}"
--         hPutStrLn stderr $ printf "%d checking, found %d" (idx :: Int) (length ss)

-- Actual functions
symmetryCandidates :: (Operator -> Bool) -> [HVec] -> Int -> Operator -> [Operator]
symmetryCandidates bptester _ 0 q 
    | bptester s = [s]
    | otherwise  = []
    where
        s = makeS q
symmetryCandidates bptester xs d q = concat $ map (symmetryCandidates bptester xs (d-1)) $ extendProjector q xs
    where
        extendProjector q xs  = map ((q +) . proj) $ xs `onlyOrthoToProj` q

makeS :: Operator -> Operator
makeS q = (ident d) - (2 `scale` q)
    where
        d = rows q

-- symmetryCandidates' :: (Operator -> Bool) -> [(HVec, Traces)] -> Int -> (Operator, Traces) -> [Operator]
-- symmetryCandidates' bptester _ 0 (q, _) 
--     | bptester s = [s]
--     | otherwise  = []
--     where
--         s = makeS q
-- symmetryCandidates' bptester xs d qq@(q, traces) = concat $ map (symmetryCandidates' bptester xs (d-1)) $ qs 
--     where
--         qs = map (extendByVector qq) $ onlyPromisingVectors xs (q, traces)
-- 
-- extendByVector :: (Operator, Traces) -> (HVec, Traces) -> (Operator, Traces)
-- extendByVector (q, t1) (x, t2) = (q + proj x, t1 + t2)
-- 
-- onlyPromisingVectors :: [(HVec, Traces)] -> (Operator, Traces) -> [(HVec, Traces)]
-- onlyPromisingVectors xs (q, traces) = filter (\(x, t) -> isValidTrace (traces + t)) $ oxs
--     where
--         oxs = onlyOrthoToProj' xs q 
--         isValidTrace t = (maximum . (map magnitude) . toList) t < k
--         k = fromIntegral (n - 1) / 2.0 
--         n = round . sqrt . fromIntegral . rows $ q
-- 
-- -- doesntBreakUnital q x = k > maxvalue
-- --     where
-- --         k = fromIntegral (n - 1) / 2.0 
-- --         n = round . sqrt . fromIntegral . rows $ q
-- --         values = map (unitalTest (q + proj x)) $ basisProjectors
-- --         maxvalue = maximum $ map magnitude values
-- --         basisProjectors = map (proj . (basis n)) $ [1..n]
-- -- 
-- --
unitalTraces :: HVec -> Traces
unitalTraces v = fromList $ map (unitalTrace $ proj v) ps
    where
        n = round . sqrt . fromIntegral . dim $ v
        ps = map (proj . (basis n))  $ [1..n]
        unitalTrace q p = trace (q <> ((ident n) `kronecker` p))

unitalTraces' :: Operator -> Traces
unitalTraces' a = fromList $ map (unitalTrace a) ps
    where
        n = round . sqrt . fromIntegral . rows $ a
        ps = map (proj . (basis n))  $ [1..n]
        unitalTrace q p = trace (q <> ((ident n) `kronecker` p))

unitalTrace q p = trace (q <> ((ident n) `kronecker` p))
    where
        n = round . sqrt . fromIntegral . rows $ q
        

-- -- unitalTest q p = trace (q <> ((ident n) `kronecker` p))
-- --     where
-- --         n = round . sqrt . fromIntegral . rows $ q

onlyOrthoToProj :: [HVec] -> Operator -> [HVec]
onlyOrthoToProj xs q = filter (isInNullSpace q) xs

onlyOrthoToProj' :: [(HVec, a)] -> Operator -> [(HVec, a)]
onlyOrthoToProj' xs q = filter (\(x, t) -> isInNullSpace q x) xs

isInNullSpace :: Operator -> HVec -> Bool
isInNullSpace q x = equalsZero $ x `cdot` (q <> x)

--
-- Combinatorics
--

-- Returns set of basic vectors with defined Schmidt rank
basicSchmidtRankSet :: Int -> Int -> [[Complex Double]] -> [HVec]
basicSchmidtRankSet rank dim phases = filter ((==rank) . schmidtRank) $ concat $ map addPairsWithPhase phases
    where
        vlist = tensorBasis dim
        sets = subsets rank vlist
        addPairsWithPhase f = map (sum . (zipWith scale f)) sets

subsets :: Int -> [a] -> [[a]]
subsets 0 xs = []
subsets 1 xs = map (\x -> [x]) xs
subsets _ [] = []
subsets n (x:xs) = (map (x:) $ subsets (n-1) xs) ++ (subsets n xs)


--
-- Tools
--

equalsZero :: Complex Double -> Bool
equalsZero n = (magnitude n) < accuracy

sumLists :: Num a => [a] -> [a] -> [a]
sumLists a b = map (\(x,y) -> x + y) $ zip a b
--
-- Formatting
--
formatMatrix :: Operator -> String
formatMatrix a = intercalate "\n" rows
       where 
             rows = map formatVector $ toRows a

formatVector :: HVec -> String
formatVector v = "[" ++ (intercalate ", " (toStrings v)) ++ "]"
       where 
             toStrings = (map formatComplex) . toList 

formatComplex :: Complex Double -> String
formatComplex z
    | (abs re < accuracy) && (abs im < accuracy) = "    0"
    | abs re < accuracy = formatReal im ++ "I"
    | abs im < accuracy = formatReal re
    | otherwise = (formatReal re) ++ " + " ++ (formatReal im) ++ "I"
    where
        re = realPart z
        im = imagPart z

formatReal :: Double -> String
formatReal x
    | abs x < accuracy = "    0"
    | isInteger x = printf "%5.0g" x
    | otherwise = printf "%5.2g" x

isInteger :: Double -> Bool
isInteger x = abs (ix - x) < accuracy
    where
        ix = fromIntegral $ round x

mfMatrix a = "{\n" ++ (intercalate ",\n" rows) ++ "}\n"
    where
        rows = map formatRow $ toRows a
        formatRow r = "{" ++ 
                        (intercalate ", " (toStrings r)) ++ 
                        "}"
        toStrings r = map mfComplex $ toList r

mfComplex z
    | (abs re < accuracy) && (abs im < accuracy) = "0"
    | abs re < accuracy = mfReal im ++ "I"
    | abs im < accuracy = mfReal re
    | otherwise = (mfReal re) ++ " + " ++ (mfReal im) ++ "I"
    where
        re = realPart z
        im = imagPart z

mfReal x
    | abs x < accuracy = "0"
    | otherwise = printf "%.2g" x
accuracy = 1e-6
imaginary = 0 :+ 1

--
-- Data
-- 
phases2C = [
           [1.0, 1.0],
           [1.0, -1.0],
           [1.0, 0.0 :+ 1.0],
           [1.0, 0.0 :+ (-1.0)]
           ] :: [[Complex Double]]
phases2R = [
           [1.0, 1.0],
           [1.0, -1.0]
           ] :: [[Complex Double]]
phases3R = [
           [1.0, 1.0, 1.0],
           [1.0, 1.0, -1.0],
           [1.0, -1.0, 1.0]
           ] :: [[Complex Double]]
c = [[1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1]] :: [[Complex Double]]
choiT3 = fromLists c

qT3 = ((ident 9) - choiT3)/2

c2 = [[-1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, -1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, -1]] :: [[Complex Double]]
choiCC = fromLists c2

c3 = [[0, 0, 0, 1],
      [0, 0, 1, 0],
      [0, 1, 0, 0],
      [1, 0, 0, 0]] :: [[Complex Double]]
mc3 = fromLists c3

(bs2:bss2) = basicSchmidtRankSet 2 2 phases2R
(bs3:bss3) = basicSchmidtRankSet 2 3 phases2R
