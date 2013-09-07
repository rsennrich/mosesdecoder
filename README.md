WARNING: this version is not being regularly updated. Use the main repository at https://github.com/moses-smt/mosesdecoder - the flexibility scores are integrated in the main repository.

This is a sample implementation of flexibility scores, as described in

Sennrich, Rico (2013): Promoting Flexible Translations in Statistical Machine Translation.
    In: Proceedings of Machine Translation Summit XIV, Nice, France.

To use, simply add the command line option --flexibility-score when calling the train-model.perl script.

The order of the scores in the phrase table will be:

p(s|t) lex(s|t) flex_left(s|t) flex_right(s|t) [flex_sub(s|t)] p(t|s) lex(t|s) flex_left(t|s) flex_right(t|s) [flex_sub(t|s)] phrase_penalty
