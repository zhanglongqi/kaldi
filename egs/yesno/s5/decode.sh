decode_cmd="utils/run.pl"
# Decoding

printf "\n*******************************Decoding *******************************\n"
steps/decode.sh --nj 1 --cmd "$decode_cmd" \
    exp/mono0a/graph_tgpr data/test_yesno exp/mono0a/decode_test_yesno
printf "\n*******************************Decoding *******************************\n"
