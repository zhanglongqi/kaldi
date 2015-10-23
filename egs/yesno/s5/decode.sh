decode_cmd="utils/run.pl"
# Decoding
. ./path.sh
#printf "\n*******************************apply-cmvn *******************************\n"

#apply-cmvn  --utt2spk=ark:data/test_yesno/split1/1/utt2spk scp:data/test_yesno/split1/1/cmvn.scp scp:data/test_yesno/split1/1/feats.scp ark:-  | add-deltas  ark:- ark:- 
#printf "\n*******************************end apply-cmvn *******************************\n"
#printf "\n*******************************add-d *******************************\n"

#add-deltas  ark:- ark:-
#printf "\n******************************end add-d *******************************\n"

printf "\n*******************************Decoding *******************************\n"

gmm-latgen-faster --max-active=7000 --beam=13.0 --lattice-beam=6.0 --acoustic-scale=0.083333 --allow-partial=true --word-symbol-table=exp/mono0a/graph_tgpr/words.txt exp/mono0a/final.mdl exp/mono0a/graph_tgpr/HCLG.fst 'ark,s,cs:apply-cmvn  --utt2spk=ark:data/test_yesno/split1/1/utt2spk scp:data/test_yesno/split1/1/cmvn.scp scp:data/test_yesno/split1/1/feats.scp ark:- | add-deltas  ark:- ark:- |' 'ark:|gzip -c > exp/mono0a/decode_test_yesno/lat.1.gz' 

printf "\n*******************************Decoding *******************************\n"
