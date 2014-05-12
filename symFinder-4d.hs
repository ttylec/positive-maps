import HPositive
import Control.Parallel.Strategies
import Control.DeepSeq
import System.Environment
import Data.List


dim = 4
bptester = randomBpTest4 1000
rank = 6
vrank = 2
phases = phases3R
aim = 3 -- (dim-1)/2 * vrank

main :: IO ()
main = do
        let
            (x:xs) = basicSchmidtRankSet vrank dim phases
        -- findSymmetries bptester rank xs
        parallelConstructSymmetry bptester aim (x:xs)
