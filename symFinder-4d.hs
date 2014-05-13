import HPositive
import Control.Parallel.Strategies
import Control.DeepSeq
import System.Environment
import Data.List


dim = 4
bptester = randomBpTest4 1000
aim = 3 -- (dim-1)/2 * vrank
rank = 3
vrank = 2
phases = phases2C

main :: IO ()
main = do
        let
            (x:xs) = basicSchmidtRankSet vrank dim phases
            q0 = proj x
            counts@(c0:_) = snd . head . addCounts $ [x]
            ss = constructSymmetry aim (q0, counts) $ addCounts xs
            bpss = filter bptester ss
    
        -- putStrLn . show $ counts
        -- putStrLn . show . length $ xs
        -- putStrLn $ show $ length $ bpss
        -- putStrLn $ show $ length $ parss
        -- parallelConstructSymmetryStepII 1 bptester aim (q0, counts) $ addCounts xs
        parallelConstructSymmetry bptester aim (x:xs)
