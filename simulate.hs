{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wall #-}
import Prelude hiding (reads, read)
import Data.Int
import Data.Random
import Data.Random.Distribution.Categorical
import Data.Random.Distribution.Uniform
import Data.Random.Source.DevRandom
import Control.Monad
import Control.Monad.IO.Class
import Options.Applicative
import Bio.Core.Sequence
import Bio.Sequence.Fasta
import qualified Data.ByteString.Lazy.Char8 as LBS
import Text.Printf

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

work :: FilePath -> FilePath -> FilePath -> Int -> Int64 -> IO ()
work infile outfile freqfile n_reads read_len =
  (flip runRVarT DevRandom :: RVarT IO () -> IO ()) $ do
    refs <- liftIO $ readFasta infile
    abundances <- simulate_abundances (length refs)
    reads <-  forM [1..n_reads] $ \i -> do
      ref :: LBS.ByteString <- categoricalT (zip abundances (map (unSD . seqdata) refs))
      start_pos <- integralUniform 0 (LBS.length ref - 1)
      let
        read0 = LBS.take read_len $ LBS.drop start_pos ref
        read = read0 <> LBS.replicate (read_len - LBS.length read0) 'A'
      return $ Seq (SeqLabel . LBS.pack $ printf "read%.4d" i) (SeqData read) Nothing
    liftIO $ do
      writeFasta outfile reads
      writeFile freqfile $ unlines $ zipWith (printf "%s\t%.5f") (map (toString . seqid) refs) abundances
