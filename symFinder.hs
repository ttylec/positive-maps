import HPositive
import Control.Parallel.Strategies
import Control.DeepSeq
import System.Environment
import Data.List


dim = 7
bptester = randomBpTest7 100
rank = 6
vrank = 3
phases = phases3R
aim = 6 -- (n-1)/2 * vrank

main :: IO ()
main = do
        let
            (x:xs) = basicSchmidtRankSet vrank dim phases
        -- findSymmetries bptester rank xs
        parallelConstructSymmetry bptester aim (x:xs)
