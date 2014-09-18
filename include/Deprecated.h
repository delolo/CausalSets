//CONNECTION:
//   /* the returned vector<bool> result has size N. inInterval[k] will be
//     *  true if i precedes k precedes j. */
//    vector<bool> inInterval(int i, int j) {
//        vector<bool> result(size);
//        int k;
//        for (k = i + 1; k < j; k++) {
//            if (causal(i, k) && causal(k, j)) {
//                result[k] = true;
//            }
//            else {
//                result[k] = false;
//            }
//        }
//        return result;
//    }
