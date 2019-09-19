import java.lang.ref.Reference;
import java.util.*;



public class Player {

    //maximizes for player X
    public int minimaxAlphaBeta(final GameState gameState, int depth, int alpha, int beta, int player) {

        Vector<GameState> possibleStates = new Vector<>();
        gameState.findPossibleMoves(possibleStates);

        int v;

        if (depth == 0 || possibleStates.size() == 0){
            return checkLines(player, gameState);
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

    /** score for whole board for the player */
    int checkLines(int thePlayer, GameState gameState) {
        int BOARD_SIZE = 4;
        int tot_score = 0;
        int result;
        int draws = 0;

        for (int row = 0; row < BOARD_SIZE; ++row)
            for (int col = 0; col < BOARD_SIZE; ++col) {
                tot_score += checkLine(thePlayer, row, col, 0, row, col, BOARD_SIZE - 1, gameState);

            }
        for (int row = 0; row < BOARD_SIZE; ++row)
            for (int layer = 0; layer < BOARD_SIZE; ++layer) {
                tot_score += checkLine(thePlayer, row, 0, layer, row, BOARD_SIZE - 1, layer, gameState);

            }
        for (int col = 0; col < BOARD_SIZE; ++col)
            for (int layer = 0; layer < BOARD_SIZE; ++layer) {
                tot_score += checkLine(thePlayer, 0, col, layer, BOARD_SIZE - 1, col, layer, gameState);

            }

        for (int row = 0; row < BOARD_SIZE; ++row) {
            tot_score += checkLine(thePlayer, row, 0, 0, row, BOARD_SIZE - 1, BOARD_SIZE - 1, gameState);

        }
        for (int col = 0; col < BOARD_SIZE; ++col) {
            tot_score += checkLine(thePlayer, 0, col, 0, BOARD_SIZE - 1, col, BOARD_SIZE - 1, gameState);

        }
        for (int layer = 0; layer < BOARD_SIZE; ++layer) {
            tot_score += checkLine(thePlayer, 0, 0, layer, BOARD_SIZE - 1, BOARD_SIZE - 1, layer, gameState);

        }

        for (int row = 0; row < BOARD_SIZE; ++row) {
            tot_score +=checkLine(thePlayer, row, 0, BOARD_SIZE - 1, row, BOARD_SIZE - 1, 0, gameState);

        }
        for (int col = 0; col < BOARD_SIZE; ++col) {
            tot_score += checkLine(thePlayer, 0, col, BOARD_SIZE - 1, BOARD_SIZE - 1, col, 0, gameState);

        }
        for (int layer = 0; layer < BOARD_SIZE; ++layer) {
            tot_score += checkLine(thePlayer, 0, BOARD_SIZE - 1, layer, BOARD_SIZE - 1, 0, layer, gameState);

        }

        tot_score +=checkLine(thePlayer, 0, 0, 0, BOARD_SIZE - 1, BOARD_SIZE - 1, BOARD_SIZE - 1, gameState);

        tot_score +=checkLine(thePlayer, 0, BOARD_SIZE - 1, 0, BOARD_SIZE - 1, 0, BOARD_SIZE - 1, gameState);

        tot_score += checkLine(thePlayer, BOARD_SIZE - 1, 0, 0, 0, BOARD_SIZE - 1, BOARD_SIZE - 1, gameState);

        tot_score += checkLine(thePlayer, BOARD_SIZE - 1, BOARD_SIZE - 1, 0, 0, 0, BOARD_SIZE - 1, gameState);

        if (thePlayer == Constants.CELL_X){
            return -1 * tot_score;
        }
        return tot_score;

    }

    public static int rowColumnLayerToCell(int row, int column, int layer)
    {
        int BOARD_SIZE = 4;
        return column + row * BOARD_SIZE + layer * BOARD_SIZE * BOARD_SIZE;
    }

    /** get score for a line */
    private int checkLine(int player, int row1, int col1, int layer1, int row2, int col2, int layer2, GameState gameState) {
        int score = 0;
        int BOARD_SIZE = 4;
        int dRow = (row2 - row1) / (BOARD_SIZE - 1);
        int dCol = (col2 - col1) / (BOARD_SIZE - 1);
        int dLayer = (layer2 - layer1) / (BOARD_SIZE - 1);
        int opponent = (player == Constants.CELL_X) ? Constants.CELL_O : Constants.CELL_X;

        int playerCells = 0, opponentCells = 0;

        for (int i = 0; i < BOARD_SIZE; ++i) {
            if (gameState.at(rowColumnLayerToCell(row1 + dRow * i, col1 + dCol * i, layer1 + dLayer * i)) == player) {
                playerCells++;
            }
            if (gameState.at(rowColumnLayerToCell(row1 + dRow * i, col1 + dCol * i, layer1 + dLayer * i)) == opponent) {
                opponentCells++;
            }
            if ((playerCells > 0) && (opponentCells > 0))
            {
                return 0;
            }
            if(opponentCells == 0)
            {
            switch (playerCells){
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
            if(playerCells == 0)
            {
                switch (opponentCells){
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
    /*
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
    }*/

    /*public int eval(final GameState gameState, int player) {

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
    }*/
    /*
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

    }*/


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
                if (current_score > best_score) {
                    best_move = nextStates.get(state_index);
                }
            } else {
                if (current_score < best_score) {
                    best_move = nextStates.get(state_index);
                }
            }
        }

        return best_move;
    }    
}
