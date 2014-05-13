import HPositive
import Control.Parallel.Strategies
import Control.DeepSeq
import System.Environment
import Data.List


dim = 5
bptester = randomBpTest5 100
aim = 4 -- (dim-1)/2 * vrank
vrank = 2
phases = phases2R

main :: IO ()
main = do
        let
            (x:xs) = basicSchmidtRankSet vrank dim phases
            q0 = proj x
            counts@(c0:_) = snd . head . addCounts $ [x]
            ss = constructSymmetry bptester aim (q0, counts) $ addCounts xs
            bpss = filter bptester ss
    
        -- putStrLn . show $ counts
        -- putStrLn . show . length $ xs
        -- putStrLn $ show $ length $ bpss
        -- findSymmetries bptester rank xs
        parallelConstructSymmetry bptester aim (x:xs)
