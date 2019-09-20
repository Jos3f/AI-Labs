import java.lang.ref.Reference;
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
        int BOARD_SIZE = 4;
        //sum of score for all lines
        int tot_score = 0;

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

        return new Score(tot_score, gameState);

    }

    //Help function taken from the GameState class
    public static int rowColumnLayerToCell(int row, int column, int layer)
    {
        int BOARD_SIZE = 4;
        return column + row * BOARD_SIZE + layer * BOARD_SIZE * BOARD_SIZE;
    }

    /** get score for a line (taken from GameState class)*/
    private int checkLine(int player, int row1, int col1, int layer1, int row2, int col2, int layer2, GameState gameState) {
        int score = 0;
        int BOARD_SIZE = 4;
        int dRow = (row2 - row1) / (BOARD_SIZE - 1);
        int dCol = (col2 - col1) / (BOARD_SIZE - 1);
        int dLayer = (layer2 - layer1) / (BOARD_SIZE - 1);
        int opponent = (player == Constants.CELL_X) ? Constants.CELL_O : Constants.CELL_X;

        int playerCells = 0, opponentCells = 0;

        //sums score for each element on the line
        for (int i = 0; i < BOARD_SIZE; ++i) {
            if (gameState.at(rowColumnLayerToCell(row1 + dRow * i, col1 + dCol * i, layer1 + dLayer * i)) == player) {
                playerCells++;
            }
            if (gameState.at(rowColumnLayerToCell(row1 + dRow * i, col1 + dCol * i, layer1 + dLayer * i)) == opponent) {
                opponentCells++;
            }
            //If no one can win here
            if ((playerCells > 0) && (opponentCells > 0))
            {
                return 0;
            }
        }

        //Scale the score depending on how many marks the player has on an "opponent-empty" line
        //scaling factor 74 (# lines on board) is used to make it preferrable to have two marks on one line rather than one mark on two separate lines
        if(opponentCells == 0)
        {
            switch (playerCells){
                case 0:
                    score += 0;
                    break;
                case 1:
                    score += 74;
                    break;
                case 2:
                    score += Math.pow(74,2);
                    break;
                case 3:
                    score += Math.pow(74,3);
                    break;
                case 4:
                    score += Math.pow(74,4);
                    break;
                default:
            }
            //Scale the score negativly depending on how many marks the opponent has on a "player-empty" line
        }else if(playerCells == 0)
        {
            switch (opponentCells){
                case 0:
                    score -= 0;
                    break;
                case 1:
                    score -= 74;
                    break;
                case 2:
                    score -= Math.pow(74,2);
                    break;
                case 3:
                    score -= Math.pow(74,3);
                    break;
                case 4:
                    score -= Math.pow(74,4);
                    break;
                default:
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

        //Determine which player's score to maximized
        int opponent = Constants.CELL_X;
        if (gameState.getNextPlayer() == Constants.CELL_X){
            opponent = Constants.CELL_O;
        }

        int depth_minus_one  = 1;
        Score best_scenario = minimaxAlphaBeta(gameState, depth_minus_one, Integer.MIN_VALUE, Integer.MAX_VALUE, gameState.getNextPlayer(), opponent,  gameState.getNextPlayer());

        return best_scenario.move;
        }



}
