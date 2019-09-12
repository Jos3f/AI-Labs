#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include "Matris.cp"
#include <vector>
#include <map>

namespace ducks
{

Matris model_1_A;
Matris model_1_B;
Matris model_1_pi;

std::vector<std::vector<int>> observation_sequences; // Bird, time
std::vector<std::vector<Matris>> model_1_alpha; // Bird, time step, state
std::vector<std::vector<double>> c_birds; //c for all birds

int time_step;
int current_round;
int TRAINING_FREQUENCY = 10;
int MIN_OBS_TRAIN = 10;

// EMovement = movement_to_shoot


Player::Player()
{
  model_1_A = Matris(5,5, {
    0.10, 0.20, 0.3, 0.15, 0.25,
    0.15, 0.10, 0.15, 0.35, 0.25,
    0.2, 0.08, 0.12, 0.25, 0.35,
    0.25, 0.35, 0.13, 0.07, 0.2,
    0.35, 0.15, 0.25, 0.12, 0.13
  });
  model_1_B = Matris(9,5, {
    0.10, 0.12, 0.21, 0.15, 0.12, 0.05, 0.08, 0.09, 0.08,
    0.15, 0.10, 0.15, 0.08, 0.15, 0.1, 0.08, 0.12, 0.07,
    0.2, 0.08, 0.12, 0.05, 0.10, 0.05, 0.15, 0.2, 0.05,
    0.15, 0.15, 0.11, 0.07, 0.2, 0.08, 0.1, 0.02, 0.12,
    0.15, 0.15, 0.14, 0.12, 0.13, 0.1, 0.05, 0.05, 0.11
  });
  model_1_pi = Matris(5,1, {
    0.12, 0.18, 0.2, 0.3, 0.2
  });

  time_step = -1;
  current_round = -1;
  // observation_sequences = vector<vector<int>>();
  // setUpRound();
}

// Do everything that is needed for a new round
void setUpRound(const GameState &pState, const Deadline &pDue){

  c_birds = std::vector<std::vector<double>>(pState.getNumBirds(), std::vector<double>());

  time_step = 0;
  observation_sequences = std::vector<std::vector<int>>();
  model_1_alpha = std::vector<std::vector<Matris>>();
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    observation_sequences.push_back(std::vector<int>());
    model_1_alpha.push_back(std::vector<Matris>());
  }
}

void train(Matris &A, Matris &B, Matris &pi, std::vector<std::vector<int>> &Os){
  return;
}

void updateAlpha(Matris &A, Matris &B, Matris &pi, std::vector<Matris> &alpha, std::vector<int> &observations, int bird_index){

  if (time_step == 0) {
    //compute alpha zero
    c_birds[bird_index].push_back(0);
    alpha.push_back(Matris(A.rows(), 1));
    for (int state = 0; state < A.rows(); state++) {
      alpha[0](state, 0) = pi(state,0)*B(observations[0],state);
      c_birds[bird_index][0] += alpha[0](state, 0);
    }
    // scale the alpha zero
    c_birds[bird_index][0] = 1.0 / c_birds[bird_index][0];
    for (int state = 0; state < A.rows(); state++) {
      alpha[0](state, 0) *= c_birds[bird_index][0];
    }

  }else{
    //compute alpha t
    alpha.push_back(Matris(A.rows(), 1));
    //    for (int t = 1; t < total_observations; t++) {
      c_birds[bird_index][time_step] = 0;
      for (int state = 0; state < A.rows(); state++) {
        alpha[time_step](state, 0) = 0;
        for (int state2 = 0; state2 < A.rows(); state2++) {
          alpha[time_step](state, 0) += alpha[time_step-1](state2, 0) * A(state,state2);
        }
        alpha[time_step](state, 0) = alpha[time_step](state, 0) * B(observations[time_step],state);
        c_birds[bird_index][time_step] += alpha[time_step](state, 0);
      }
      // scale the alpha 1 to T
      c_birds[bird_index][time_step] = 1.0 / c_birds[bird_index][time_step];

      for (int state = 0; state < A.rows(); state++) {
        alpha[time_step](state, 0) *= c_birds[bird_index][time_step];
      }
  }
  //Matris test_alpha = Matris(5,1, {
    //0.12, 0.18, 0.2, 0.3, 0.2
  //});
  //alpha.push_back(test_alpha);
  return;
}

//update at each time_step
void updateObservations(const GameState &pState, std::vector<std::vector<int>> &observation_sequences){
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    Bird bird = pState.getBird((int) bird_index);

    if (bird.isDead()) {
      return;
    }
    observation_sequences[(int) bird_index].push_back(bird.getLastObservation());
  }
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */
     std::cerr << "here 1" << '\n';
     time_step++;
     std::cerr << "Round: " << pState.getRound() << '\n';

     if (current_round != pState.getRound()) {
       // Do everything that is needed for a new round
       current_round = pState.getRound();
       setUpRound(pState, pDue);
     }
     updateObservations(pState, observation_sequences);
     // std::cerr << "Model 1 A" << '\n';
     // std::cerr << model_1_A << '\n';
     // std::cerr << "Model 1 B" << '\n';
     // std::cerr << model_1_B << '\n';
     // std::cerr << "Model 1 pi" << '\n';
     // std::cerr << model_1_pi << '\n';
     // std::cerr << pState.getNumBirds() << '\n';
     if (time_step % TRAINING_FREQUENCY == 0 && time_step > MIN_OBS_TRAIN) {
       train(model_1_A, model_1_B, model_1_pi, observation_sequences);
     }

     //update alpha
     std::vector<Matris> pi_birds(pState.getNumBirds(), Matris(0,0));
     std::vector<Matris> movement_predictions_birds(pState.getNumBirds(), Matris(0,0));
     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       if (pState.getBird((int) bird_index).isAlive()) {
         updateAlpha(model_1_A, model_1_B, model_1_pi,
           model_1_alpha[(int) bird_index],
           observation_sequences[(int) bird_index],
         (int) bird_index);
        pi_birds[(int) bird_index] = model_1_alpha[(int) bird_index][time_step];
        movement_predictions_birds[(int) bird_index] = pi_birds[(int) bird_index] * model_1_A * model_1_B;
          }
       }

     std::map<int, EMovement> movement_map;
     movement_map[0] = MOVE_UP_LEFT;
     movement_map[1] = MOVE_UP;
     movement_map[2] = MOVE_UP_RIGHT;
     movement_map[3] = MOVE_LEFT;
     movement_map[4] = MOVE_STOPPED;
     movement_map[5] = MOVE_RIGHT;
     movement_map[6] = MOVE_DOWN_LEFT;
     movement_map[7] = MOVE_DOWN;
     movement_map[8] = MOVE_DOWN_RIGHT;
     movement_map[9] = COUNT_MOVE;

     //make decision
     int most_certain_bird = -1;
     double threshold = 0.1;
     double probability = threshold;
     int movement_to_shoot = -1;
     for (int bird = 0; bird < pState.getNumBirds(); bird++) {
       if (pState.getBird((int) bird).isAlive()) {
         for (int movement = 0; (int) movement < (movement_predictions_birds[(int) bird]).cols(); movement++) {
           double movement_prob = (movement_predictions_birds[(int) bird])(movement, 0);

          //std::cerr << "prob" << probability<< '\n';
           if (probability < movement_prob) {
             probability = movement_prob;
             most_certain_bird = (int) bird;
             movement_to_shoot = movement;
           }
         }
       }

     }

     if (most_certain_bird != -1) {
       std::cerr << "Time step: " << time_step << '\n';
       std::cerr << "bird to shoot " << most_certain_bird << '\n';
       std::cerr << "birds pred move " << movement_map[movement_to_shoot] << '\n';
       std::cerr << "is bird to shoot alive" << pState.getBird(most_certain_bird).isAlive() <<'\n';
       return Action(most_certain_bird, movement_map[movement_to_shoot]);
     }


    // This line choose not to shoot
    return cDontShoot;

    //This line would predict that bird 0 will move right and shoot at it
    // return Action(0, MOVE_RIGHT);
}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */

    std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);
    return lGuesses;
}

void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */
    std::cerr << "HIT BIRD!!!" << std::endl;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */
}


} /*namespace ducks*/
