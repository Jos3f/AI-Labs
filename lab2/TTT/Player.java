import java.util.*;

public class Player {

    //maximizes for player X
    public int minimaxAlphaBeta(final GameState gameState, int depth, int alpha, int beta, int player) {

        Vector<GameState> possibleStates = new Vector<>();
        gameState.findPossibleMoves(possibleStates);


        int v;

        if (depth == 0 || possibleStates.size() == 0){
            return heuristic(gameState, player);
        } else if (player == Constants.CELL_X){
            v = Integer.MIN_VALUE;
            for (GameState child: possibleStates
                 ) {
                int v_child = minimaxAlphaBeta(child, depth - 1, alpha, beta, Constants.CELL_O);
                if (v_child > v){
                    v = v_child;
                }
                if (v > alpha){
                    alpha = v;
                }
                if (beta <= alpha){
                    return v;
                }
            }
        } else {
            v = Integer.MAX_VALUE;
            for (GameState child: possibleStates
            ) {
                int v_child = minimaxAlphaBeta(child, depth - 1, alpha, beta, Constants.CELL_X);
                if (v_child < v){
                    v = v_child;
                }
                if (v < beta){
                    beta = v;
                }
                if (beta <= alpha){
                    return v;
                }
            }
        }
        return v;
    }


    // Maybe change name
    public int heuristic(final GameState gameState, int player) {
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

        if (player == Constants.CELL_X){
            return -1 * score;
        }
        return score;
    }

    public int overallMarks(int [] lines, int [] lines_opponent){
        int score = 0;

        for (int line = 0; line < lines.length; line++) {
            if (lines_opponent[line] == 0){
                switch (lines[line]){
                    case 0:
                        score += 0;
                        break;
                    case 1:
                        score += 1;
                        break;
                    case 2:
                        score += 5;
                        break;
                    case 3:
                        score += 25;
                        break;
                    case 4:
                        score += 1250;
                        break;
                    default:
                }
            }
            if (lines[line] == 0){
                switch (lines_opponent[line]){
                    case 0:
                        score -= 0;
                        break;
                    case 1:
                        score -= 1;
                        break;
                    case 2:
                        score -= 5;
                        break;
                    case 3:
                        score -= 25;
                        break;
                    case 4:
                        score -= 1250;
                        break;
                    default:
                }
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

        int depth = 4;
        GameState best_move = nextStates.get(0);
        int best_score = minimaxAlphaBeta(nextStates.get(0), depth, Integer.MIN_VALUE, Integer.MAX_VALUE, (nextStates.get(0)).getNextPlayer());


        for (int state_index = 1; state_index < nextStates.size(); state_index++) {
            int current_score = minimaxAlphaBeta(nextStates.get(state_index), depth, Integer.MIN_VALUE, Integer.MAX_VALUE, (nextStates.get(state_index)).getNextPlayer());

            if (gameState.getNextPlayer() == Constants.CELL_X) {
                if (current_score < best_score) {
                    best_move = nextStates.get(state_index);
                }
            } else {
                if (current_score > best_score) {
                    best_move = nextStates.get(state_index);
                }
            }
        }

        return best_move;
        /**
         * Here you should write your algorithms to get the best next move, i.e.
         * the best next state. This skeleton returns a random move instead.
         */
        // Random random = new Random();
        // return nextStates.elementAt(random.nextInt(nextStates.size()));
    }    
}
