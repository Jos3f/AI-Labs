import java.util.*;

public class Player {

    //class to hold moves and its heuristic function value
    private class Score {

        public int score_val;
        public GameState move;

        public Score(int score_val, GameState move){
            this.score_val = score_val;
            this.move = move;
        }

    }

    //maximizes for maxI_player
    public Score minimaxAlphaBeta(final GameState gameState, int depth, int alpha, int beta, int player, int opponent, int maxi_player) {

        Vector<GameState> possibleStates = new Vector<>();
        gameState.findPossibleMoves(possibleStates);


        Score v = new Score(0, null);

        //Base case
        if (depth == 0 || possibleStates.size() == 0){
            return heuristic(maxi_player, gameState);


        } else if (player == maxi_player){ //Recursive case for maxi_player
            v.score_val = Integer.MIN_VALUE;
            for (GameState child: possibleStates
            ) { //   nent and player swith at next level in tree
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


    // Maybe change name
    public Score heuristic( int player, final GameState gameState) {
        int [] rows = {0,0,0,0};
        int [] cols = {0,0,0,0};
        int [] diag = {0,0};

        int [] rows_opponent = {0,0,0,0};
        int [] cols_opponent = {0,0,0,0};
        int [] diag_opponent = {0,0};

        int opponent = Constants.CELL_O;
        if (player == Constants.CELL_O){
            opponent = Constants.CELL_X;
        }

        for (int row = 0; row < rows.length ; row++) {
            for (int col = 0; col < cols.length; col++) {
                if (gameState.at(row, col) == player) {
                    rows[row]++;
                    cols[col]++;
                    if (row == col) { // Diagonal right down
                        diag[0]++;
                    } else if ((row + col) == 3) { // Diagonal left down
                        diag[1]++;
                    }
                } else if (gameState.at(row, col) == opponent) {
                    rows_opponent[row]++;
                        cols_opponent[col]++;
                        if (row == col) { // Diagonal right down
                            diag_opponent[0]++;
                        } else if ((row + col) == 3) { // Diagonal left down
                            diag_opponent[1]++;
                        }
                }
            }
        }

        int score = 0;

        // Calc myMarks for rows
        int row_score = overallMarks(rows, rows_opponent);
        // Calc myMarks for cols
        int col_score = overallMarks(cols, cols_opponent);
        // Calc myMarks for diag
        int diag_score = overallMarks(diag, diag_opponent);

        score = row_score + col_score + diag_score;

        Score return_object = new Score(score, gameState);
   /*
        if (player == Constants.CELL_X){
            return score;
        }
        return -1 * score;
        */
        return return_object;
    }

    public int overallMarks(int [] lines, int [] lines_opponent){
        int score = 0;

        for (int line = 0; line < lines.length; line++) {
            if (lines_opponent[line] == 0 && lines[line] == 0){
                continue;
            } else if (lines_opponent[line] == 0){
                score += Math.pow(5, lines[line]);
            } else if (lines[line] == 0){
                score -= Math.pow(5, lines_opponent[line]);
            }
        }
        return score;
    }

    /**
     * Performs a move
     *
     * @param gameState
     *            the current state of the board
     * @param deadline
     *            time before which we must have returned
     * @return the next state the board is in after our move
     */
    public GameState play(final GameState gameState, final Deadline deadline) {
        Vector<GameState> nextStates = new Vector<GameState>();
        gameState.findPossibleMoves(nextStates);


        if (nextStates.size() == 0) {
            // Must play "pass" move if there are no other moves possible.
            return new GameState(gameState, new Move());
        }

        int opponent = Constants.CELL_O;
        if (gameState.getNextPlayer() == Constants.CELL_O){
            opponent = Constants.CELL_X;
        }
        int depth = 5;

        Score best_suggestion = minimaxAlphaBeta(gameState, depth, Integer.MIN_VALUE, Integer.MAX_VALUE, gameState.getNextPlayer(), opponent, gameState.getNextPlayer());

        return best_suggestion.move;
        /**
         * Here you should write your algorithms to get the best next move, i.e.
         * the best next state. This skeleton returns a random move instead.
         */
        // Random random = new Random();
        // return nextStates.elementAt(random.nextInt(nextStates.size()));
    }    
}
