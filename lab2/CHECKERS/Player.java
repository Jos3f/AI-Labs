import java.util.*;

public class Player {

    public State_dictionary global_dictionary;

    private class Pieces_state {
        public int red_pawns;
        public int red_kings;
        public int white_pawns;
        public int white_kings;

        // Initialize pieces state to default board
        public Pieces_state() {
            red_pawns = 12;
            red_kings = 0;
            white_pawns = 12;
            white_kings = 0;
        }

        // Copy object
        public Pieces_state(Pieces_state p) {
            red_pawns = p.red_pawns;
            red_kings = p.red_kings;
            white_pawns = p.white_pawns;
            white_kings = p.white_kings;
        }
    }

    /**
     * Used for repeated state checking
     */
    private class State_keeper {
        public int[] board;
        public int next_player;
        public  int[][] table;
        public int depth;

        public State_keeper(GameState gameState, int depth , int[][] zobrist_init_table){
            this.next_player = gameState.getNextPlayer();
            board = new int[gameState.NUMBER_OF_SQUARES];
            for (int square_index = 0; square_index < gameState.NUMBER_OF_SQUARES; square_index++) {
                board[square_index] = gameState.get(square_index);
            }
            table = zobrist_init_table;
            this.depth = depth;
        }

        @Override
        public int hashCode() {
            int hash = 0;
            for (int board_index = 0; board_index < board.length; board_index++) {
                if ((board[board_index]) != 0){ // Zero means empty cell
                    int piece = board[board_index];
                    // System.err.println(piece);
                    hash = hash ^ table[board_index][piece];
                }
            }
            if (next_player == 2){
                hash = hash * -1;
            }
            hash += depth;
            return hash;
        }

        @Override
        public boolean equals(Object obj) {
            return true;
        }
    }

    private class State_dictionary {
        HashMap<State_keeper, Integer> dictionary;
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

        public int getScore(GameState gameState, int depth){
            State_keeper state = new State_keeper(gameState, depth, zobrist_table);
            return dictionary.get(state);
        }

        public void setScore(GameState gameState, int depth,  int score){
            State_keeper state = new State_keeper(gameState, depth, zobrist_table);
            dictionary.put(state, score);
        }

        public boolean stateExists(GameState gameState, int depth){
            State_keeper state = new State_keeper(gameState, depth, zobrist_table);
            return dictionary.get(state) != null;
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
            return heuristicFunction(maxi_player, gameState);


        } else if (player == maxi_player){ //Recursive case for maxi_player
            v.score_val = Integer.MIN_VALUE;
            for (GameState child: possibleStates
            ) { //opponent and player switch at next level in tree
                Score v_child;
                if (! global_dictionary.stateExists(child, depth)){
                    v_child = minimaxAlphaBeta(child, depth - 1, alpha, beta, opponent, player, maxi_player);
                    global_dictionary.setScore(child, depth, v_child.score_val);
                } else {
                    v_child = new Score(0, null);
                    v_child.score_val = global_dictionary.getScore(child, depth);
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
                    v_child.score_val = global_dictionary.getScore(child, depth);
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
         *   value of my pawn on my side    = 3;
         *   value of my pawn on their side = 6;
         *   value of my king               = 9;
         * */
        int opponent = Constants.CELL_RED;
        if (gameState.getNextPlayer() == Constants.CELL_RED){
            opponent = Constants.CELL_WHITE;
        }
/*
        int my_pawns = getNumberOf(gameState, thePlayer, 0);
        int my_kings = getNumberOf(gameState, thePlayer, Constants.CELL_KING);
        int opponent_pawns = getNumberOf(gameState, opponent, 0);
        int opponent_kings = getNumberOf(gameState, opponent, Constants.CELL_KING);
    */

        int[] value_of_cell = new int[]{
                1, 1, 1, 1,
                3, 2, 2, 2,
                3, 3, 3, 4,
                5, 4, 4, 4,
                5, 5, 5, 6,
                7, 6, 6, 6,
                7, 7, 7, 8,
                9, 9, 9, 9,
        };

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




        for (int pos = 0; pos < GameState.NUMBER_OF_SQUARES; pos++) {
            if ((gameState.get(pos) & thePlayer) > 0){ // My piece
                if ((gameState.get(pos) & Constants.CELL_KING) > 0){ // My king piece
                    score += 8 * value_of_cell_king[pos];
                } else { // My pawn piece
                    score += 4 * value_of_cell[pos];

                }

            } else if ((gameState.get(pos) & opponent) > 0){ // opponent piece
                if ((gameState.get(pos) & Constants.CELL_KING) > 0){ // opponent king piece
                    score -= 8 * value_of_cell_king[(GameState.NUMBER_OF_SQUARES - 1) - pos];

                } else { // opponent pawn piece
                    score -= 4 * value_of_cell[(GameState.NUMBER_OF_SQUARES - 1) - pos];

                }
            }
        }
/*
        score+= 4 * my_pawns;
        score+= 8 * my_kings;
        score-= 4 * opponent_pawns;
        score-= 8 * opponent_kings;
*/


        return new Score(score, gameState);
    }

    /**
     * Function for calculating average distances between the pieces
     */
    public double avgDistance(GameState gameState){
        ArrayList<Integer> my_pieces = new ArrayList<>(); // Fill with pos of my pieces
        ArrayList<Integer> opponet_pieces = new ArrayList<>(); // Fill with pos of opponent pieces
        double dist = 0;
         return dist;
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

/*
        // Testing hashfunction
        dict.setScore(pState, 10,1337);
        boolean stateExist = dict.stateExists(pState, 10);
        System.err.println("State exists: " + stateExist);
        stateExist = dict.stateExists(lNextStates.get(0), 10);
        System.err.println("State exists: " + stateExist);

        int[] a = new int[]{1,2};
        int[] b = new int[]{1,2};
        HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
        map.put(a.hashCode(), 123);
        System.err.println("Key b gives: " + map.get(b.hashCode()));

        System.exit(1);
*/

        Score best_scenario = new Score(0,null);
        int depth_minus_one = 1;

        while (pDue.timeUntil() > 750000000) {
            global_dictionary = new State_dictionary();

            best_scenario = minimaxAlphaBeta(pState, depth_minus_one, Integer.MIN_VALUE, Integer.MAX_VALUE, pState.getNextPlayer(), opponent, pState.getNextPlayer());

            depth_minus_one++;
        }

        /**
         * Here you should write your algorithms to get the best next move, i.e.
         * the best next state. This skeleton returns a random move instead.
         */
        System.err.println("NUMBER OF RED PAWNS: " + getNumberOf(pState, Constants.CELL_RED, 0));
        System.err.println("NUMBER OF RED KINGS: " + getNumberOf(pState, Constants.CELL_RED, Constants.CELL_KING));
        System.err.println("NUMBER OF WHITE PAWNS: " + getNumberOf(pState, Constants.CELL_WHITE, 0));
        System.err.println("NUMBER OF WHITE KINGS: " + getNumberOf(pState, Constants.CELL_WHITE, Constants.CELL_KING));

        return (best_scenario.move);
    }
}
