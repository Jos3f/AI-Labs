#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include "Matris.cp"
#include <vector>
#include <map>

namespace ducks
{

// Forward declarations
void updateAlpha(Matris &A, Matris &B, Matris &pi, std::vector<Matris> & alpha, std::vector<int> &observations, std::vector<double> & c);
void updateAlphaTrain(Matris &A, Matris &B, Matris &pi, std::vector<Matris> & alpha, std::vector<int> &observations, std::vector<double> & c);

Matris model_1_A;
Matris model_1_B;
Matris model_1_pi;

std::vector<Matris> model_bird_A;
std::vector<Matris> model_bird_B;
std::vector<Matris> model_bird_pi;


std::vector<std::vector<int>> observation_sequences; // Bird, time
std::vector<std::vector<Matris>> model_1_alpha; // Bird, time step, state
std::vector<std::vector<Matris>> model_1_beta; // Bird, time step, state
std::vector<Matris> model_1_alpha_train; // time step, state
std::vector<double> c_train; // time step, state
std::vector<std::vector<double>> c_birds; //c for all birds,

int time_step;
int current_round;
int TRAINING_FREQUENCY = 10;
int MIN_OBS_TRAIN = 10;

// EMovement = movement_to_shoot

struct HMM {
  Matris a;
  Matris b;
  Matris pi;

  // int num_states;
  // int num_observations;

  int total_observations;

  // double [] c;

  HMM(){
    a = Matris(5,5, {
      0.10, 0.20, 0.3, 0.15, 0.25,
      0.15, 0.10, 0.15, 0.35, 0.25,
      0.2, 0.08, 0.12, 0.25, 0.35,
      0.25, 0.35, 0.13, 0.07, 0.2,
      0.35, 0.15, 0.25, 0.12, 0.13
    });

    b = Matris(9,5, {
      0.10, 0.12, 0.21, 0.15, 0.12, 0.05, 0.08, 0.09, 0.08,
      0.15, 0.10, 0.15, 0.08, 0.15, 0.1, 0.08, 0.12, 0.07,
      0.2, 0.08, 0.12, 0.05, 0.10, 0.05, 0.15, 0.2, 0.05,
      0.15, 0.15, 0.11, 0.07, 0.2, 0.08, 0.1, 0.02, 0.12,
      0.15, 0.15, 0.14, 0.12, 0.13, 0.1, 0.05, 0.05, 0.11
    });

    pi = Matris(5,1, {
      0.12, 0.18, 0.2, 0.3, 0.2
    });

    // num_states = A.rows();
    // num_observations = B.cols();

  }

  Matris calcNextEmissionPD(Matris current_pi){
    return current_pi * a * b;
  }

  // std::vector<Matris> calcAlpha(std::vector<int> &observations){
  //
  //   return;
  // }

  void train(std::vector<int> observations){

  }

};

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
  model_1_beta = std::vector<std::vector<Matris>>();
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    observation_sequences.push_back(std::vector<int>());
    model_1_alpha.push_back(std::vector<Matris>());
    model_1_beta.push_back(std::vector<Matris>());
  }
}

void train(Matris &A, Matris &B, Matris &pi, std::vector<Matris> alpha, std::vector<int> & observation_sequences, std::vector<double> & c){
  // updateAlphaTrain(A, B, pi,
  //   alpha,
  //   observation_sequences, c);
  // updateBeta(model_1_A, model_1_B,
  //   model_1_beta[(int) bird_index],
  //   observation_sequences[(int) bird_index],
  //   c_birds[(int) bird_index]);
  return;
}

void updateAlpha(Matris &A, Matris &B, Matris &pi, std::vector<Matris> & alpha, std::vector<int> &observations, std::vector<double> & c){
  int observation_time_step = 0;
  if (time_step == 0) {
    //compute alpha zero
    c.push_back(0);
    alpha.push_back(Matris(A.rows(), 1));
    for (int state = 0; state < A.rows(); state++) {
      alpha[0](state, 0) = pi(state,0)*B(observations[0],state);
      c[0] += alpha[0](state, 0);
    }
    // scale the alpha zero
    c[0] = 1.0 / c[0];
    for (int state = 0; state < A.rows(); state++) {
      alpha[0](state, 0) *= c[0];
    }
  } else {
    //compute alpha t
    alpha.push_back(Matris(A.rows(), 1));
    c.push_back(0);
      c[time_step] = 0;
      for (int state = 0; state < A.rows(); state++) {
        (alpha[time_step])(state, 0) = 0;
        for (int state2 = 0; state2 < A.rows(); state2++) {
          (alpha[time_step])(state, 0) += alpha[time_step-1](state2, 0) * A(state,state2);
        }
        (alpha[time_step])(state, 0) = alpha[time_step](state, 0) * B(observations[time_step],state);
        c[time_step] += alpha[time_step](state, 0);
      }
      // scale the alpha 1 to T
      c[time_step] = 1.0 / c[time_step];
      for (int state = 0; state < A.rows(); state++) {
        (alpha[time_step])(state, 0) *= c[time_step];
      }
  }
  return;
}

void updateAlphaTrain(Matris &A, Matris &B, Matris &pi, std::vector<Matris> & alpha, std::vector<int> &observations, std::vector<double> & c){
  int observation_time_step = observations.size();
  // Reset alpha and c
  alpha = std::vector<Matris>();
  c = std::vector<double>();
  std::cerr << "Alpha: " << '\n';

  //compute alpha zero
  c.push_back(0);
  alpha.push_back(Matris(A.rows(), 1));
  for (int state = 0; state < A.rows(); state++) {
    alpha[0](state, 0) = pi(state,0)*B(observations[0],state);
    c[0] += alpha[0](state, 0);
  }
  // scale the alpha zero
  c[0] = 1.0 / c[0];
  for (int state = 0; state < A.rows(); state++) {
    alpha[0](state, 0) *= c[0];
  }

  //compute alpha t
  for (int t = 1; t < observation_time_step; t++) {

    alpha.push_back(Matris(A.rows(), 1));
    c.push_back(0);
    c[time_step] = 0;
    for (int state = 0; state < A.rows(); state++) {
      (alpha[t])(state, 0) = 0;
      for (int state2 = 0; state2 < A.rows(); state2++) {
        (alpha[t])(state, 0) += alpha[t-1](state2, 0) * A(state,state2);
      }
      (alpha[t])(state, 0) = alpha[t](state, 0) * B(observations[t],state);
      c[t] += alpha[t](state, 0);
    }
    // scale the alpha 1 to T
    c[t] = 1.0 / c[t];
    for (int state = 0; state < A.rows(); state++) {
      (alpha[t])(state, 0) *= c[t];
    }
  }

  return;
}

void updateBeta(Matris & A, Matris & B,  std::vector<Matris> & beta, std::vector<int> &observations, std::vector<double> & c){
  // calculate Beta pass T
  beta.insert(beta.begin(), Matris(A.rows(),1));
  for (int state = 0; state < A.rows(); state++) {
    beta[0](state, 0) = 1 * c[observations.size() - 1];
  }

  for (int observation_step = observations.size() - 2; observation_step >= 0; observation_step--) {
    beta.insert(beta.begin(), Matris(A.rows(),1));
    for (int state = 0; state < A.rows(); state++) {
      for (int state2 = 0; state2 < A.rows(); state2++) {
        beta[0](state, 0) += A(state2, state) * B(observations[observation_step + 1], state2) * beta[1](state2, 0);
      }
      // Scale beta t with same factor as alpha t
      beta[0](state, 0) *= c[observation_step];
      // allBeta(observation_step, i) *= c[observation_step];
    }
  }

}

//update at each time_step
void updateObservations(const GameState &pState, std::vector<std::vector<int>> &observation_sequences){
  // std::cerr << "update observation call" << '\n';
  // std::cerr << "Num of birds: " << pState.getNumBirds() <<  '\n';
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    Bird bird = pState.getBird((int) bird_index);

    // std::cerr << "For bird index: " << bird_index << " :" << bird.getLastObservation() << '\n';
    if (bird.isAlive()) {
      observation_sequences[(int) bird_index].push_back(bird.getLastObservation());
    }
  }
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */
     // std::cerr << "here 1" << '\n';
     time_step++;
     std::cerr << "Time step start: " << time_step << '\n';

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
       std::vector<int> observation_sequences_train = std::vector<int>();
       for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
         observation_sequences_train.insert(observation_sequences_train.end(), observation_sequences[(int) bird_index].begin(), observation_sequences[(int) bird_index].begin());
       }
       model_1_alpha_train = std::vector<Matris>();
       c_train = std::vector<double>();
       train(model_1_A, model_1_B, model_1_pi, model_1_alpha_train, observation_sequences_train, c_train);
     }

     //update alpha and beta
     std::vector<Matris> pi_birds(pState.getNumBirds(), Matris(0,0));
     std::vector<Matris> movement_predictions_birds(pState.getNumBirds(), Matris(0,0));
     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       if (pState.getBird((int) bird_index).isAlive()) {
         updateAlpha(model_1_A, model_1_B, model_1_pi,
           model_1_alpha[(int) bird_index],
           observation_sequences[(int) bird_index],
           c_birds[(int) bird_index]);

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

     // most_certain_bird = 5;
     if (most_certain_bird != -1) {
       std::cerr << "Time step: " << time_step << '\n';
       std::cerr << "bird to shoot " << most_certain_bird << '\n';
       std::cerr << "birds pred move " << movement_map[movement_to_shoot] << '\n';
       std::cerr << "Current score " << pState.myScore() <<'\n';
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
