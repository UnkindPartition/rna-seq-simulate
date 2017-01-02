{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wall #-}
import Prelude hiding (reads, read)
import Data.Int
import Data.Maybe
import Data.Random
import Data.Random.Distribution.Categorical
import Data.Random.Distribution.Uniform
import Data.Random.Source.DevRandom
import qualified Data.Map as Map
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.State
import Options.Applicative
import Bio.Core.Sequence
import Bio.Sequence.Fasta
import qualified Data.ByteString.Lazy.Char8 as LBS
import Text.Printf

-- | Return random relative (unnormalized) number of transcripts
simulate_abundances :: Int -> RVarT m [Double]
simulate_abundances n = do
  replicateM n $ exp <$> normalT (-2) 1

seq_length :: Sequence -> Double
seq_length = fromIntegral . LBS.length . unSD . seqdata

calc_observed_tpm :: Map.Map SeqLabel Double -> [Sequence] -> [Double]
calc_observed_tpm freq_map seqs =
  let
    rpk = map (\s -> (fromMaybe 0 $ Map.lookup (seqid s) freq_map) / seq_length s) seqs
  in map (* (1e6 / sum rpk)) rpk

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
    transcript_abundances <- simulate_abundances (length refs)
    let nucleotide_abundances = zipWith (\f r -> f * seq_length r) transcript_abundances refs
    -- A Map in the State monad contain the actual frequencies of reads
    -- from each isoform. These frequencies are generated from the
    -- categorical distribution based on abundances.
    (reads, freqs) <-  flip runStateT Map.empty $
      forM [1..n_reads] $ \i -> do
        ref <- lift $ weightedCategoricalT (zip nucleotide_abundances refs)
        modify' $ Map.insertWith (+) (seqid ref) 1
        let ref_seq = (unSD . seqdata) ref
        start_pos <- lift $ integralUniform 0 (LBS.length ref_seq - 1)
        let
          read0 = LBS.take read_len $ LBS.drop start_pos ref_seq
          read = read0 <> LBS.replicate (read_len - LBS.length read0) 'A'
        return $ Seq (SeqLabel . LBS.pack $ printf "read%.4d" i) (SeqData read) Nothing
    let
      true_tpm = map (* (1e6  / sum transcript_abundances)) transcript_abundances
      observed_tpm = calc_observed_tpm freqs refs
    liftIO $ do
      writeFasta outfile reads
      writeFile freqfile $ unlines $
        "seq_id\ttrue_tpm\tobserved_tpm" :
        zipWith3 (printf "%s\t%.4f\t%.4f") (map (toString . seqid) refs)
          true_tpm
          observed_tpm
