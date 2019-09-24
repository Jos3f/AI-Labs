import java.util.*;

public class Player {

    //GLobal dictionary for storing visited states (to be able to evaluate more possible states)
    public State_dictionary global_dictionary;

    /**
     * Used for repeated state checking
     * Hashes board using zobrist to a key to use in hash map
     */
    private class State_keeper {
        public int[] board;
        public int next_player;
        public  int[][] table;

        public State_keeper(GameState gameState, int[][] zobrist_init_table){
            this.next_player = gameState.getNextPlayer();
            board = new int[gameState.NUMBER_OF_SQUARES];
            //inserts content of each square in board-vector
            for (int square_index = 0; square_index < gameState.NUMBER_OF_SQUARES; square_index++) {
                board[square_index] = gameState.get(square_index);
            }
            table = zobrist_init_table;
        }

        //board -> key using Zobrist hashing
        @Override
        public int hashCode() {
            int hash = 0;
            for (int board_index = 0; board_index < board.length; board_index++) {
                if ((board[board_index]) != 0){ // Zero means empty cell
                    int piece = board[board_index];
                    hash = hash ^ table[board_index][piece];
                }
            }
            //to make unique hashes for each player. Matters who's turn it is
            if (next_player == 2){
                hash = hash * -1;
            }
            return hash;
        }
        // check if state_keeper objects are equal. Hashmap uses this to confirm that it is the same object
        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;

            if (obj == null || obj.getClass() != this.getClass())
                return false;

            State_keeper state = (State_keeper) obj;

            if (this.next_player != state.next_player)
                return false;
            //checks if boards are the same
            for (int pos = 0; pos < board.length; pos++) {
                if (this.board[pos] != state.board[pos]){
                    return false;
                }
            }
            return true;
        }
    }

    // class for keeping all the states
    private class State_dictionary {
        HashMap<State_keeper, Score_and_depth> dictionary;
        int[][] zobrist_table; // Zobris initial map

        public State_dictionary() {
            dictionary = new HashMap<>();

            // Initialize Zobrist table
            zobrist_table = new int[GameState.NUMBER_OF_SQUARES][7];
            Random rand = new Random();
            for (int cell_index = 0; cell_index < zobrist_table.length; cell_index++) {
                for (int piece_type = 1; piece_type < (zobrist_table[0]).length; piece_type++) {
                    zobrist_table[cell_index][piece_type] = rand.nextInt();
                }
            }
        }

        //gets score for a state
        public Score_and_depth getScore(GameState gameState){
            State_keeper state = new State_keeper(gameState, zobrist_table);
            return dictionary.get(state);
        }

        // sets score abd depth as values
        public void setScore(GameState gameState, int depth,  int score){
            State_keeper state = new State_keeper(gameState, zobrist_table);
            dictionary.put(state, new Score_and_depth(score, depth));
        }

        // if stored depth is equal to or greater than current deptband boards are same;
        // we wanna use the stored value (=> return state already exists)
        public boolean stateExists(GameState gameState, int depth){
            State_keeper state = new State_keeper(gameState, zobrist_table);
            Score_and_depth old_score = dictionary.get(state);
            if (old_score == null){
                return false;
            }
            if (old_score.depth >= depth){
                return true;
            }
            return false;
        }
    }

    //class to hold moves and its heuristic function value
    private class Score {

        public int score_val;
        public GameState move;

        public Score(int score_val, GameState move){
            this.score_val = score_val;
            this.move = move;
        }

    }

    public int getNumberOf(GameState pState, int color, int pawn_or_king){
        int count = 0;
        int looking_for = color | pawn_or_king;

        for (int pos = 0; pos < GameState.NUMBER_OF_SQUARES; pos++) {
            if (pState.get(pos) == looking_for){
                count++;
            }
        }
        return count;
    }

    //maximizes for maxI_player
    public Score minimaxAlphaBeta(final GameState gameState, int depth, int alpha, int beta, int player, int opponent, int maxi_player) {

        Vector<GameState> possibleStates = new Vector<>();
        gameState.findPossibleMoves(possibleStates);


        Score v = new Score(0, null);

        //Base case
        if (depth == 0 || possibleStates.size() == 0){
            //calculate value always from the same player's perspective
            return heuristicFunction(maxi_player, gameState);


        } else if (player == maxi_player){ //Recursive case for maxi_player
            v.score_val = Integer.MIN_VALUE;
            for (GameState child: possibleStates
            ) { //opponent and player switch at next level in tree
                Score v_child;
                if (! global_dictionary.stateExists(child, depth)){ //compute and add/overwrite in hashmap
                    v_child = minimaxAlphaBeta(child, depth - 1, alpha, beta, opponent, player, maxi_player);
                    global_dictionary.setScore(child, depth, v_child.score_val);
                } else {
                    v_child = new Score(0, null);
                    v_child.score_val = global_dictionary.getScore(child).score;
                }
                if (v_child.score_val > v.score_val){
                    v.score_val = v_child.score_val;
                    v.move = child;
                }
                if (v.score_val > alpha){
                    alpha = v.score_val;
                }
                if (beta <= alpha){
                    return v;
                }
            }
        } else { //Recursive case for non maxi_player
            v.score_val = Integer.MAX_VALUE;
            for (GameState child: possibleStates
            ) {
                Score v_child;
                if (! global_dictionary.stateExists(child, depth)){
                    v_child = minimaxAlphaBeta(child, depth - 1, alpha, beta, opponent, player, maxi_player);
                    global_dictionary.setScore(child, depth, v_child.score_val);
                } else {
                    v_child = new Score(0, null);
                    v_child.score_val = global_dictionary.getScore(child).score;
                }
                if (v_child.score_val < v.score_val){
                    v.score_val = v_child.score_val;
                    v.move = child;
                }
                if (v.score_val < beta){
                    beta = v.score_val;
                }
                if (beta <= alpha){
                    return v;
                }
            }
        }
        return v;
    }

    /** computes value and move for whole board for the player (obtained form the GameState class)*/
    Score heuristicFunction(int thePlayer, GameState gameState) {
        int score = 0;
        /**
         *
         * Calculate the score in this block
         *   basic value of my pawn    = 1;
         *   basic value of my king    = 9;
         * */
        int opponent = Constants.CELL_RED;
        if (gameState.getNextPlayer() == Constants.CELL_RED){
            opponent = Constants.CELL_WHITE;
        }

        // Special case when game is over (only cares about who won or draw)
        if (gameState.isEOG()){
            if (gameState.isDraw()){
                return new Score(0, gameState);
            } else if (gameState.isRedWin()){
                if (thePlayer == Constants.CELL_RED){
                    return new Score(999999999, gameState);
                } else {
                    return new Score(-999999999, gameState);
                }
            } else if (gameState.isWhiteWin()){
                if (thePlayer == Constants.CELL_WHITE){
                    return new Score(999999999, gameState);
                } else {
                    return new Score(-999999999, gameState);
                }
            }
        }

        // better to have pawns near opponents boarder and on the sides
        int[] value_of_cell = new int[]{
                1, 1, 1, 1,
                3, 2, 2, 2,
                3, 3, 3, 4,
                5, 4, 4, 4,
                5, 5, 5, 6,
                7, 6, 6, 6,
                7, 7, 7, 8,
                9, 9, 9, 9,  // pawn will not be here since it will be a king by then
        };

        // kings are equally valuable all over the board vertically, since they can move backwards
        int[] value_of_cell_king = new int[]{
                2, 2, 2, 2,
                2, 1, 1, 1,
                1, 1, 1, 2,
                2, 1, 1, 1,
                1, 1, 1, 2,
                2, 1, 1, 1,
                1, 1, 1, 2,
                2, 2, 2, 2,
        };

        // iterate though every cell and assign score based on heuristics
        for (int pos = 0; pos < GameState.NUMBER_OF_SQUARES; pos++) {

            if ((gameState.get(pos) & thePlayer) > 0){ // My piece
                if ((gameState.get(pos) & Constants.CELL_KING) > 0){ // My king piece
                    score += 9 * value_of_cell_king[pos];
                } else { // My pawn piece
                    score += 1 * value_of_cell[pos];

                }

            } else if ((gameState.get(pos) & opponent) > 0){ // opponent piece
                if ((gameState.get(pos) & Constants.CELL_KING) > 0){ // opponent king piece
                    score -= 9 * value_of_cell_king[(GameState.NUMBER_OF_SQUARES - 1) - pos];

                } else { // opponent pawn piece
                    score -= 1 * value_of_cell[(GameState.NUMBER_OF_SQUARES - 1) - pos];

                }
            }
        }
        return new Score(score, gameState);
    }

    /**
     * Performs a move
     *
     * @param pState
     *            the current state of the board
     * @param pDue
     *            time before which we must have returned
     * @return the next state the board is in after our move
     */
    public GameState play(final GameState pState, final Deadline pDue) {

        Vector<GameState> lNextStates = new Vector<GameState>();
        pState.findPossibleMoves(lNextStates);

        if (lNextStates.size() == 0) {
            // Must play "pass" move if there are no other moves possible.
            return new GameState(pState, new Move());
        }

        int opponent = Constants.CELL_RED;
        if (pState.getNextPlayer() == Constants.CELL_RED){
            opponent = Constants.CELL_WHITE;
        }

        State_dictionary dict = new State_dictionary();

        /**
         * Unit test 2020 (copyright) below
         */
/*
        // Testing hashfunction
        dict.setScore(pState, 10,1337);
        GameState pState_copy = pState;
        boolean stateExist = dict.stateExists(pState_copy,10);
        System.err.println("State exists: " + stateExist);
        stateExist = dict.stateExists(lNextStates.get(0), 10);
        System.err.println("State exists: " + stateExist);

       // int[] a = new int[]{1,2};
     //   int[] b = new int[]{1,2};
   //     HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
 //       map.put(a, 123);
//        System.err.println("Key b gives: " + map.get(b));

        System.exit(1);
*/

        Score best_scenario = new Score(0,null);
        int depth_minus_one = 1;


        while (pDue.timeUntil() > 600000000) {
            global_dictionary = new State_dictionary();

            if (best_scenario.score_val >= 999999999){
                return best_scenario.move;
            }

            best_scenario = minimaxAlphaBeta(pState, depth_minus_one, Integer.MIN_VALUE, Integer.MAX_VALUE, pState.getNextPlayer(), opponent, pState.getNextPlayer());

            depth_minus_one++;
        }

        //Debugging
        System.err.println("NUMBER OF RED PAWNS: " + getNumberOf(pState, Constants.CELL_RED, 0));
        System.err.println("NUMBER OF RED KINGS: " + getNumberOf(pState, Constants.CELL_RED, Constants.CELL_KING));
        System.err.println("NUMBER OF WHITE PAWNS: " + getNumberOf(pState, Constants.CELL_WHITE, 0));
        System.err.println("NUMBER OF WHITE KINGS: " + getNumberOf(pState, Constants.CELL_WHITE, Constants.CELL_KING));

        return (best_scenario.move);
    }

    // Class to hold score and depth used as value in hash map
    private class Score_and_depth {
        public int score;
        public int depth;

        Score_and_depth(int s, int d){
            score = s;
            depth = d;
        }
    }
}

