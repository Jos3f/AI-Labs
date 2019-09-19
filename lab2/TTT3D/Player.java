import java.lang.ref.Reference;
import java.util.*;




public class Player {

    //maximizes for player X
    public int minimaxAlphaBeta(final GameState gameState, int depth, int alpha, int beta, int player) {

        Vector<GameState> possibleStates = new Vector<>();
        gameState.findPossibleMoves(possibleStates);

        /*for (GameState state: possibleStates
        ) {

            System.err.println(state.toString(state.getNextPlayer()));
            System.err.println("------------------------------");
        }*/

        /*
        int opponent = Constants.CELL_O;
        if (player == Constants.CELL_O){
            opponent = Constants.CELL_X;
        }
*/
        int v;

        if (depth == 0 || possibleStates.size() == 0){
            return eval(gameState, player);
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

    public int eval(final GameState gameState, int player) {

        int opponent = Constants.CELL_O;
        if (player == Constants.CELL_O){
            opponent = Constants.CELL_X;
        }

        ArrayList<int []> game_rows_player = insertGameRows(gameState, player);
        ArrayList<int []> game_rows_opponent = insertGameRows(gameState, opponent);

        int score = 0;

        // Calc myMarks for rows
        for (int line_group_index = 0; line_group_index < game_rows_player.size(); line_group_index++) {
            score += overallMarks(game_rows_player.get(line_group_index),
                    game_rows_opponent.get(line_group_index));
        }

        if (player == Constants.CELL_X){
            return -1 * score;
        }
        return score;
    }

    public ArrayList<int []>  insertGameRows(GameState gameState, int player){
        int [] rows1 =      {0,0,0,0,
                0,0,0,0,
                0,0,0,0,
                0,0,0,0};

        int [] rows2 =      {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
        int [] rows3 =      {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
        int [] diag_row1 =  {0,0,0,0, 0,0,0,0};
        int [] diag_row2 =  {0,0,0,0, 0,0,0,0};
        int [] diag_row3 =  {0,0,0,0, 0,0,0,0};
        int [] diag_main =  {0,0,0,0};


        for (int layer1 = 0; layer1 < 4; layer1++) {
            for (int layer2 = 0; layer2 < 4; layer2++) {
                for (int layer3 = 0; layer3 < 4; layer3++) {
                    if (gameState.at(layer1, layer2, layer3) == player) {

                        rows1[layer1 + 4 * layer2]++;
                        rows2[layer2 + 4 * layer3]++;
                        rows3[layer3 + 4 * layer1]++;

                        // Diagonals for rows1
                        if (layer1 == layer2) { // Diagonal right down
                            diag_row1[layer3]++;
                        } else if ((layer1 + layer2) == 3) { // Diagonal left down
                            diag_row1[layer3 + 4]++;

                        }
                        // Diagonals for rows2
                        if (layer2 == layer3) { // Diagonal right down
                            diag_row2[layer1]++;
                        } else if ((layer2 + layer3) == 3) { // Diagonal left down
                            diag_row2[layer1 + 4]++;
                        }
                        // Diagonals for rows3
                        if (layer3 == layer1) { // Diagonal right down
                            diag_row3[layer2]++;
                        } else if ((layer3 + layer1) == 3) { // Diagonal left down
                            diag_row3[layer2 + 4]++;
                        }

                        // Main diagonals
                        if ((layer1 == layer2) && (layer1 == layer3)) {
                            diag_main[0]++;
                        } else if (layer1 == layer2 && (layer1 + layer3 == 3)) {
                            diag_main[1]++;
                        } else if ((layer1 + layer2) == 3 && (layer1 == layer3)) {
                            diag_main[2]++;
                        } else if ((layer1 + layer2) == 3 && (layer1 + layer3 == 3)) {
                            diag_main[3]++;
                        }

                    }
                }
            }
        }
        return new ArrayList<int []> (List.of(rows1, rows2, rows3, diag_row1, diag_row2, diag_row3, diag_main));

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
                        score -= 5;
                        break;
                    case 2:
                        score -= 25;
                        break;
                    case 3:
                        score -= 125;
                        break;
                    case 4:
                        score -= 6250;
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

        int depth = 1;
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
    }    
}
