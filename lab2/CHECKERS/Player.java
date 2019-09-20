import java.util.*;

public class Player {

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
            ) { //opponent and player swith at next level in tree
                Score v_child = minimaxAlphaBeta(child, depth - 1, alpha, beta, opponent, player, maxi_player);
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
                Score v_child = minimaxAlphaBeta(child, depth - 1, alpha, beta, opponent, player, maxi_player);
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

        int my_pawns = getNumberOf(gameState, thePlayer, 0);
        int my_kings = getNumberOf(gameState, thePlayer, Constants.CELL_KING);
        int opponent_pawns = getNumberOf(gameState, opponent, 0);
        int opponent_kings = getNumberOf(gameState, opponent, Constants.CELL_KING);

        score+= 4 * my_pawns;
        score+= 8 * my_kings;
        score-= 4 * opponent_pawns;
        score-= 8 * opponent_kings;

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

        int depth_minus_one = 6;
        Score best_scenario = minimaxAlphaBeta(pState, depth_minus_one, Integer.MIN_VALUE, Integer.MAX_VALUE, pState.getNextPlayer(), opponent, pState.getNextPlayer());

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