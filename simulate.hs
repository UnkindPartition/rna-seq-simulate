{-# OPTIONS_GHC -Wall #-}
import Prelude hiding (reads, read)
import Data.Random
import Data.Random.Distribution.Categorical
import Data.Random.Distribution.Uniform
import Data.Random.Source.DevRandom
import Control.Monad
import Control.Monad.IO.Class
import Options.Applicative

simulate_abundances :: Int -> RVarT m [Double]
simulate_abundances n = do
  xs <- replicateM n $ exp <$> normalT (-2) 1
  return $ map (/ sum xs) xs

main :: IO ()
main = join . execParser $
  info (helper <*> parser)
  (  fullDesc
  <> header "A simple RNA-Seq read simulator"
  <> progDesc "Generate reads based on a reference transcriptome"
  )
  where
    parser :: Parser (IO ())
    parser =
      work
        <$> strOption   (long "infile"   <> short 'i' <> metavar "FILE")
        <*> strOption   (long "outfile"  <> short 'o' <> metavar "FILE")
        <*> strOption   (long "freqfile" <> short 'f' <> metavar "FILE")
        <*> option auto (long "nreads"   <> short 'n' <> metavar "NUMBER")
        <*> option auto (long "readlen"  <> short 'l' <> metavar "NUMBER")

work :: FilePath -> FilePath -> FilePath -> Int -> Int -> IO ()
work infile outfile freqfile n_reads read_len =
  (flip runRVarT DevRandom :: RVarT IO () -> IO ()) $ do
    refs <- lines <$> liftIO (readFile infile)
    abundances <- simulate_abundances (length refs)
    reads <-  replicateM n_reads $ do
      ref <- categoricalT $ zip abundances refs
      start_pos <- integralUniform 0 (length ref - 1)
      let read = take read_len $ drop start_pos ref ++ repeat 'A' 
      return read
    liftIO $ do
      writeFile outfile $ unlines reads
      writeFile freqfile $ unlines reads
