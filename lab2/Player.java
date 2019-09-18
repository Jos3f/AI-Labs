import java.util.*;

public class Player {

    public int minimaxAlphaBeta(final GameState gameState, int depth, int alpha, int beta, int player) {

        Vector<GameState> possibleStates = new Vector<>();
        gameState.findPossibleMoves(possibleStates);


        if (depth == 0 || possibleStates.size() == 0){
            return eval(gameState, player);
        } else if (player == Constants.CELL_X){
            int v = 0;
            for (GameState state: possibleStates
                 ) {

            }
        }
        return minimaxAlphaBeta(gameState, depth, alpha, beta, player);
    }

    public int eval(final GameState gameState, int player) {
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
                    if (gameState.at(row, col) == player) {
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
        }

        int score = 0;

        // Calc myMarks for rows
        int row_score = myMarks(rows);
        // Calc myMarks for cols
        int col_score = myMarks(cols);
        // Calc myMarks for diag
        int diag_score = myMarks(diag);

        score = row_score + col_score + diag_score;

        return score;
    }

    public int myMarks(int [] lines){
        int score = 0;
        for (int line = 0; line < lines.length; line++) {
            score += lines[line];
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

        /**
         * Here you should write your algorithms to get the best next move, i.e.
         * the best next state. This skeleton returns a random move instead.
         */
        Random random = new Random();
        return nextStates.elementAt(random.nextInt(nextStates.size()));
    }    
}
