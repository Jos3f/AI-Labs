#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include "Matris.cp"
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <time.h>
// #pragma GCC optimize ("O0")
// #include <r>


namespace ducks
{

// Forward declarations
struct HMM;
std::vector<double> guessSpecies(std::vector<std::vector<double>> lastAlphas);
int getMin(const std::vector<double> & v);
int getMax(const std::vector<double> & v);
void scale(std::vector<double> & v);
void updateObservations(const GameState &pState, std::vector<std::vector<int>> &observation_sequences);

// int guessSpeciesAndBlackProb(std::vector<int> observations, double & black_stork_prob);


// Global mappings
std::map<int, EMovement> movement_map;
std::map<int, ESpecies> species_map;

std::vector<std::vector<int>> observation_sequences; // Bird, time
std::vector<std::vector<int>> observation_sequences_species; // species, time
std::vector<int> species_observed; // Counter for how many times a species has occured so far
int birds_observed;
std::vector<double> prob_species_prior;

std::vector<std::vector<std::vector<double>>> allAlphas; // birds in round, species, pattern/hidden
std::vector<std::vector<double>> prob_species;          // Birds in round, species

int time_step;
int current_round;
int count_bird_hits = 0;
int count_shooting_attempts = 0;

int MAX_ITER_TRAIN = 100;
int TRAINING_START = 70;
double SHOOT_PROB_THRESHOLD = 0.35;
double SHOOT_STORK_PROB_THRESHOLD = 0.25;
int NUM_STATES = 5;

std::vector<HMM> AllModels;


struct HMM {
  Matris a;
  Matris b;
  Matris pi;

  int num_states;
  // int num_observations;

  int total_observations;

  std::vector<double> c;

  HMM(int num_states){

    a = generate_a_matris(num_states);/*Matris(5,5, {
      0.10, 0.20, 0.3, 0.15, 0.25,
      0.15, 0.10, 0.15, 0.35, 0.25,
      0.2, 0.08, 0.12, 0.25, 0.35,
      0.25, 0.35, 0.13, 0.07, 0.2,
      0.35, 0.15, 0.25, 0.12, 0.13
    });*/



    b = generate_b_matris(num_states);/*Matris(9,5, {
      0.10, 0.12, 0.21, 0.15, 0.12, 0.05, 0.08, 0.09, 0.08,
      0.15, 0.10, 0.15, 0.08, 0.15, 0.1, 0.08, 0.12, 0.07,
      0.2, 0.08, 0.12, 0.05, 0.10, 0.05, 0.15, 0.2, 0.05,
      0.15, 0.15, 0.11, 0.07, 0.2, 0.08, 0.1, 0.02, 0.12,
      0.15, 0.15, 0.14, 0.12, 0.13, 0.1, 0.05, 0.05, 0.11
    });*/



    pi = generate_pi(num_states);/*Matris(5,1, {
      0.12, 0.18, 0.2, 0.3, 0.2
    });*/


    // num_states = A.rows();
    // num_observations = B.cols();

  }

  double generate_element_val(int num_cols){

      double uni_rand = (double)rand() / RAND_MAX;
      double scale_factor =  (double) num_cols;
      double avg_element_value = 1/scale_factor;
      double min_val = (avg_element_value - avg_element_value/2);
      double max_val = (avg_element_value + avg_element_value/2);

      double element_val = min_val + uni_rand * (max_val - min_val);
      return element_val;
  }

  Matris generate_matris(int num_cols, int num_states){
      Matris matris(num_cols,num_states);
      for (int row = 0; row < num_states; ++row) {
          double row_sum = 0.0;

          for (int col = 0; col < num_cols; ++col) {
              double generated_element_val = generate_element_val(num_cols);
              row_sum += generated_element_val;
              matris(col,row) = generated_element_val;

          }
          for (int col = 0; col < num_cols; ++col) {
              matris(col,row) = matris(col,row)/row_sum;
          }
      }
      return matris;

  }

  Matris generate_a_matris(int num_states){
      a = generate_matris(num_states, num_states);

      return a;
  }

  Matris generate_b_matris(int num_states){
      b = generate_matris(9, num_states);

      return b;
  }
  Matris generate_pi(int num_states){
      pi = generate_matris(num_states, 1);

      return pi;
  }

  std::vector<double> calcNextAlpha(){

  }

  Matris calcNextObsDist(Matris current_pi){
    return current_pi * a * b;
  }

  Matris calcAlpha(const std::vector<int> &observations){
    int total_observations = observations.size();
    Matris allAlpha(total_observations, a.rows());
    c = std::vector<double>(total_observations, 0); // This has been changed scince HMM3
    c[0] = 0;
    for (int i = 0; i < a.rows(); i++) {
      allAlpha(0, i) = pi(i,0)*b(observations[0],i);
      c[0] += allAlpha(0, i);
    }
    // scale the alpha zero
    c[0] = 1.0 / c[0];

    // TEMPORARY OBS
    if (isinf(c[0]) || isnan(c[0])) {
      c[0] = 10e+50;
    }

    for (int i = 0; i < a.rows(); i++) {
      allAlpha(0, i) *= c[0];
    }

    //compute alpha t
    for (int t = 1; t < total_observations; t++) {
      c[t] = 0;
      for (int i = 0; i < a.rows(); i++) {
        allAlpha(t, i) = 0;
        for (int j = 0; j < a.rows(); j++) {
          allAlpha(t, i) += allAlpha(t-1, j) * a(i,j);
        }
        allAlpha(t, i) = allAlpha(t, i) * b(observations[t],i);
        c[t] += allAlpha(t, i);
      }

      // scale the alpha 1 to T
      c[t] = 1.0 / c[t];

      // TEMPORARY OBS
      if (isinf(c[t]) || isnan(c[t]) ) {
        c[t] = 10e+50;
      }

      for (int i = 0; i < a.rows(); i++) {
        allAlpha(t, i) *= c[t];
      }
    }

    return allAlpha;
  }

  Matris calcBeta(const std::vector<int> &observations){
    int total_observations = observations.size();
    Matris allBeta(total_observations, a.rows());
    // calculate Beta pass T
    for (int state = 0; state < a.rows(); state++) {
      allBeta(total_observations - 1, state) = 1 * c[total_observations - 1];
    }

    // std::cout << "Allbeta: " << '\n';
    // std::cout << allBeta << '\n';

    for (int time_step = total_observations - 2; time_step >= 0; time_step--) {
      for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.rows(); j++) {
          allBeta(time_step, i) += a(j, i) * b(observations[time_step + 1], j) * allBeta(time_step + 1, j);
        }
        // std::cout << "Beta: " << allBeta(time_step, i) << '\n';
        // Scale beta t with same factor as alpha t
        allBeta(time_step, i) *= c[time_step];
      }
    }
    return allBeta;
  }

  void smoothRow(int row, Matris & matris){
     double max_val = 0.0;
     int max_col = -1;
     for (int col = 0; col < matris.cols(); col++) {
         if (matris(col, row) > max_val){
            max_val = matris(col, row);
            max_col = col;
        }
       }
     for (int col = 0; col < matris.cols(); col++) {
         if (matris(col, row) < 0.0001){
           matris(col, row) += 0.0001;
           matris(max_col, row) -= 0.0001;
       }
     }
   }

  void train(const std::vector<int> & observations){
    double oldLogProb = -std::numeric_limits<double>::infinity();
    double logProb = -std::numeric_limits<double>::infinity();
    int maxIter = MAX_ITER_TRAIN;
    int itertation = 0;

    while ((logProb > oldLogProb || itertation == 0) && itertation < maxIter) {

      Matris allAlpha = calcAlpha(observations);
      Matris allBeta = calcBeta(observations);

      int total_observations = observations.size();

      //Compute sum of alpha T
      double denom = sumAlphaT(allAlpha);
      double resultProb = denom;
      //compute gamma i,j and gamma i
      Matris digamma[total_observations];
      Matris gamma(total_observations,a.rows());

      for (int time_step = 0; time_step < total_observations - 1; time_step++) {
        digamma[time_step] = Matris(a.rows(), a.rows());
        for (int i = 0; i < a.rows(); i++) {
          for (int j = 0; j < a.rows(); j++) {
            digamma[time_step](j,i) = allAlpha(time_step, i)*a(j,i)*b(observations[time_step+1], j)*allBeta(time_step+1, j)/denom;
            gamma(time_step, i) += digamma[time_step](j,i);
          }
        }
      }

      // special case for gamma T
      denom = 0;
      for (int i = 0; i < a.rows(); i++) {
        denom += allAlpha(total_observations-1, i);
      }
      for (int i = 0; i < a.rows(); i++) {
        gamma(total_observations-1, i) = allAlpha(total_observations-1, i)/denom;
      }

      // Re-estimate A, B and pi

      // Re-estimate pi
      for (int i = 0; i < a.rows(); i++) {
        pi(i,0) = gamma(0, i);
      }

      // Re-estimate A
      for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.rows(); j++) {
          double numer = 0;
          double denom = 0;
          for (int time_step = 0; time_step < total_observations - 1; time_step++) {
            numer = numer + digamma[time_step](j, i);
            denom = denom + gamma(time_step, i);
          }
          // std::cout << "Estimate A step("<< i << " " << j <<  "): denom: " << denom << " numer: " << numer << '\n';
          // if (denom < 0.0000000000001 && denom > -0.0000000000001 ) {
          //  a(j,i) = 0;
          // } else {

            a(j,i) = numer/denom;
          // }
        }
      }

      // Re-estimate B
      for (int j = 0; j < a.rows(); j++) {
        for (int k = 0; k < b.cols(); k++) {
          double numer = 0;
          double denom = 0;
          for (int time_step = 0; time_step < total_observations - 1; time_step++) {
            if (observations[time_step] == k) {
              numer = numer + gamma(time_step, j);
            }
            denom = denom + gamma(time_step, j);
          }
          // if (denom < 0.000000000000001 && denom > -0.0000000000001 ) {
          //   b(k,j) = 0;
          // } else {
            b(k,j) = numer/denom;
          // }
        }
      }

      // Compute log[P(O|λ)]
      oldLogProb = logProb;
      logProb = 0;
      for (int time_step = 0; time_step < total_observations; time_step++) {
        logProb -= log(c[time_step]);
      }
      itertation++;
    }

    // for (int x = 0; x < a.cols(); x++) {
    //   for (int y = 0; y < a.rows(); y++) {
    //     if (a(x,y) < 0.0001) {
    //       a(x,y) = 0.0001;
    //     }
    //   }
    // }

    for (int row = 0; row < a.rows(); row++) {
      smoothRow(row, a);
    }

    // for (int x = 0; x < b.cols(); x++) {
    //   for (int y = 0; y < b.rows(); y++) {
    //     if (b(x,y) < 0.0001) {
    //       b(x,y) = 0.0001;
    //     }
    //   }
    // }

    for (int row = 0; row < b.rows(); row++) {
      smoothRow(row, b);
    }

    // for (int x = 0; x < pi.cols(); x++) {
    //   for (int y = 0; y < pi.rows(); y++) {
    //     if (pi(x,y) < 0.0001) {
    //       pi(x,y) = 0.0001;
    //     }
    //   }
    // }

    for (int row = 0; row < pi.rows(); row++) {
      smoothRow(row, pi);
    }

    return;
  }

  Matris getLastAlpha(const std::vector<int> & observations){
    Matris allAlpha = calcAlpha(observations);
    Matris lastAlpha(a.rows() ,1);
    for (int i = 0; i < a.rows(); i++) {
      lastAlpha(i, 0) = allAlpha(observations.size() - 1, i);
    }
    return lastAlpha;
  }

  double sumAlphaT(const Matris allAlpha){
    double sum = 0;
    for (int state = 0; state < a.rows(); state++)
      sum += allAlpha(allAlpha.cols() - 1, state);
    return sum;
  }

  double calcProb (const std::vector<int> observations){
    Matris allAlpha = calcAlpha(observations);
    // std::cerr << "A matrix: " << '\n';
    // std::cerr << a << '\n';
    // std::cerr << "B matrix: " << '\n';
    // std::cerr << b << '\n';
    // std::cerr << "pi matrix: " << '\n';
    // std::cerr << pi << '\n';
    // std::cerr << "allAlpha: " << allAlpha << "\n";
    // std::cerr << '\n';
    // std::cerr << "C: ";
    // for (auto  i: c) {
    //   std::cerr << i << ", ";
    // }
    // std::cerr << '\n';
    int total_observations = observations.size();
    double logProb;
    for (int time_step = 0; time_step < total_observations; time_step++) {
      logProb -= log(c[time_step]);
    }

    return logProb;
  }

};

Player::Player()
{
  // srand(time(NULL));

  current_round = -1;
  species_observed = std::vector<int>(COUNT_SPECIES, 0);
  prob_species_prior = std::vector<double>(COUNT_SPECIES, 0);
  // observation_sequences = vector<vector<int>>();
  // setUpRound();
  AllModels = std::vector<HMM>();
  for (int i = 0; i < COUNT_SPECIES; i++) {
    AllModels.push_back(HMM(NUM_STATES));
  }

  observation_sequences_species = std::vector<std::vector<int>>();
  for (size_t species = 0; species < COUNT_SPECIES; species++) {
    observation_sequences_species.push_back(std::vector<int>());
  }

  // Init global mapping
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

  species_map[-1] = SPECIES_UNKNOWN;
  species_map[0] = SPECIES_PIGEON;
  species_map[1] = SPECIES_RAVEN;
  species_map[2] = SPECIES_SKYLARK;
  species_map[3] = SPECIES_SWALLOW;
  species_map[4] = SPECIES_SNIPE;
  species_map[5] = SPECIES_BLACK_STORK;

}

// Do everything that is needed for a new round
void setUpRound(const GameState &pState, const Deadline &pDue){

  time_step = -1;
  prob_species = std::vector<std::vector<double>>(pState.getNumBirds(), std::vector<double>(COUNT_SPECIES, 0));

  allAlphas = std::vector<std::vector<std::vector<double>>>(pState.getNumBirds(), std::vector<std::vector<double>>(COUNT_SPECIES, std::vector<double>(NUM_STATES, 0)));
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    for (size_t species = 0; species < COUNT_SPECIES; species++) {
      for (size_t state = 0; state < NUM_STATES; state++) {
        allAlphas[bird_index][species][state] = AllModels[species].b(pState.getBird(bird_index).getLastObservation(), state) * AllModels[species].pi(state, 0);

      }


    }



    prob_species[bird_index] = guessSpecies(allAlphas[bird_index]);

  }


  // missed_birds = std::vector<int>(pState.getNumBirds());

  observation_sequences = std::vector<std::vector<int>>();
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    observation_sequences.push_back(std::vector<int>());
  }
  updateObservations(pState, observation_sequences);
}

void calcAlphas(GameState pState){
  // for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
  //   for (size_t species = 0; species < COUNT_SPECIES; species++) {
  //     for (size_t state = 0; state < NUM_STATES; state++) {
  //       allAlphas[bird_index][species][state] = AllModels[species].b(pState.getBird(bird_index).getLastObservation(), state) * AllModels[species].pi(state, 0);
  //
  //     }
  //   }
  //   prob_species[bird_index] = guessSpecies(bird_index);
  // }
}

void calcSpeciesProbAndAlphas(GameState pState){
  calcAlphas(pState);
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

Matris getMovementDistForBirdModel(int num_states, std::vector<int> observations){
    HMM birdModel = HMM(num_states);
    birdModel.train(observations);
    Matris lastAlpha = birdModel.getLastAlpha(observations);
    Matris movementDist = birdModel.calcNextObsDist(lastAlpha);

    return movementDist;
}

void addMovementDistToAvg(Matris movementDist, Matris & avgMovementDist, int num_state_models){
    for (int move = 0; move < avgMovementDist.cols(); ++move) {
        avgMovementDist(move, 0) += movementDist(move, 0)/num_state_models;
    }
}
Matris computeAvgMovementDist(std::vector<int> observations){
    //create multiple HMM:s for each bird and then take avg prediction
    std::vector<double> states_bird_model = {3,4,5};
    Matris avgMovementDist(COUNT_MOVE, 1);
    int num_models = states_bird_model.size();

    //avg movement dist from stored models
    /*
    if (time_step > 90 && current_round > 7) {
        num_models += AllModels[likely_species].size();
        for (int species_model_idx = 0; species_model_idx < AllModels[likely_species].size(); ++species_model_idx) {
            HMM hmm_model = AllModels[likely_species][species_model_idx];
            Matris lastAlpha = hmm_model.getLastAlpha(observations);
            Matris movementDist = hmm_model.calcNextObsDist(lastAlpha);
            addMovementDistToAvg(movementDist, avgMovementDist, num_models);
        }
    }*/

    //Compute movement dist for each nr of states for a bird
    for (int index = 0; index < states_bird_model.size(); ++index) {
        Matris movementDist = getMovementDistForBirdModel(states_bird_model[index], observations);
        addMovementDistToAvg(movementDist, avgMovementDist, num_models);
    }
    //std::cerr << "avgMovementDist " << avgMovementDist << '\n';
    return avgMovementDist;
}
//
// bool isProbablyBlackStork(const std::vector<double> & guess_distribution){
//
//   int most_certain_bird = getMin(guess_distribution);
//
//   if (guess_distribution[most_certain_bird] / guess_distribution[SPECIES_BLACK_STORK] > SHOOTING_TOP_SIMILARITY_THRESHOLD) {
//     return true;
//   }
//   return false;
//   /*
//   double top_a = 1;
//   double top_b = 1;
//   int top_a_index = -1;
//   int top_b_index = -1;
//
//   for (size_t probability_index = 0; probability_index < guess_distribution.size(); probability_index++) {
//     if ((guess_distribution[probability_index] < top_a) || (guess_distribution[probability_index] < top_b)) {
//       if (top_a > top_b) {
//         top_a = guess_distribution[probability_index];
//         top_a_index = probability_index;
//       }else{
//         top_b = guess_distribution[probability_index];
//         top_b_index = probability_index;
//       }
//     }
//   }
//
//   std::cerr << "GUESS DIST: " << '\n';
//   for (auto& val : guess_distribution) {
//     std::cerr << val << ", ";
//   }
//   std::cerr  << '\n';
//
//
//
//   if (top_a_index == SPECIES_BLACK_STORK || top_b_index == SPECIES_BLACK_STORK) {
//     return true;
//   }
//   else {
//     return false;
//   }*/
//
//
// }
//
//
// bool speciesIsUnKnown(int best_guess_index, std::vector<double> & guess_distribution){
//   //if the species is not seen before, it is unknown
//   return (AllModels[best_guess_index].size() == 0 || guess_distribution[best_guess_index] > SPECIES_UNKNOWN_THRESHOLD);
// }
//
// double getLowestLogProb()
//
// std::vector<int> getDeathPenalty(const GameState & pState){
//   int predMove = -1;
//   int guiltyBird = -1;
//   double highestProb = -1;
//
//   for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
//     if (pState.getBird(bird_index).isAlive()) {
//       std::vector<int> observations = observation_sequences[(int) bird_index];
//       //double black_stork_prob = -100000;
//       // Get log, use it as threshold when we don't have the stork in All models
//
//
//       std::vector<double> guess_distribution = guessSpecies(observation_sequences[(int) bird_index]);
//
//       if (AllModels[SPECIES_BLACK_STORK].size() == 0 ) {
//         int most_probable_species = getMax(guess_distribution);
//         if (guess_distribution[most_probable_species] < GUESS_LOGPROB_THRESHOLD) {
//           continue;
//         }
//       }
//         scale(guess_distribution);
//         int best_guess_index = getMin(guess_distribution);
//         bool species_is_unknown = speciesIsUnKnown(best_guess_index, guess_distribution);
//         bool might_be_stork = isProbablyBlackStork(guess_distribution);
//
//
//       // int likely_species = best_guess_index;
//
//       //conditions for death penalty
//       // if ((likely_species != SPECIES_BLACK_STORK)
//       // && (black_stork_prob < BLACK_STORK_LOGPROB_THRESHOLD)
//       // && (likely_species != SPECIES_UNKNOWN) ) {
//
//       if ( !might_be_stork && (!species_is_unknown || true)
//       && (missed_birds[bird_index] < MAX_BULLETS_PER_BIRD)
//       && (missed_species[best_guess_index] < MAX_BULLETS_PER_SPECIES || current_round < 5)) {
//         /* code */
//
//           Matris avgMovementDist = computeAvgMovementDist(observations);
//
//           //std::cerr << "avgMovementDist: " << avgMovementDist << '\n';
//
//
//           for (int movement = 0; movement < avgMovementDist.cols(); movement++) {
//           if (avgMovementDist(movement, 0) > highestProb && avgMovementDist(movement, 0) > PROB_THRESHOLD) {
//             highestProb = avgMovementDist(movement, 0);
//             guiltyBird = bird_index;
//             predMove = movement;
//           }
//         }
//       }
//     }
//   }
//
//   return {guiltyBird, predMove};
// }

std::vector<double> alphasNextStep(HMM model, std::vector<double> alphas, int observation){
  std::vector<double> alphas_next_step(NUM_STATES, 0);
  for (size_t state_from = 0; state_from < NUM_STATES; state_from++) {
    for (size_t state_to = 0; state_to < NUM_STATES; state_to++) {
      alphas_next_step[state_to] += alphas[state_from] * model.a(state_to, state_from)*model.b(observation, state_to);
    }
  }
  return alphas_next_step;
}

void nextAlphaUpdate(GameState pState){
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    if (pState.getBird(bird_index).isAlive()) {
      for (size_t species = 0; species < COUNT_SPECIES; species++) {
        allAlphas[bird_index][species] = alphasNextStep(AllModels[species], allAlphas[bird_index][species], observation_sequences[bird_index][observation_sequences[bird_index].size() - 1]);
        // std::cerr << "All alphas: ";
        // for (auto prob1:  allAlphas[bird_index][species]) {
        //   std::cerr << prob1 << ", ";
        // }
        // std::cerr << '\n';
      }
    }
  }
  // exit(0);
}

void setUpTimeStep(GameState pState){
  nextAlphaUpdate(pState);

  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    prob_species[bird_index] = guessSpecies(allAlphas[bird_index]);
  }
  return;
}

double observationProb(std::vector<double> lastAlphaVector){
  double prob = 0;
  for (size_t state = 0; state < NUM_STATES; state++) {
    prob += lastAlphaVector[state];
  }
  return prob;
}

double getMoveProb(int bird_index, int movement){
  double probability = 0.000;
  for (size_t species = 0; species < COUNT_SPECIES; species++) {
    double prob_next = observationProb(alphasNextStep(AllModels[species], allAlphas[bird_index][species], movement));
    double prob_this_time_step = observationProb(allAlphas[bird_index][species]);
    if (prob_this_time_step > 0) {
      probability += ((prob_next / prob_this_time_step) * prob_species[bird_index][species]);
    }
  }
  return probability;
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */


     if (current_round != pState.getRound()) {
       // Do everything that is needed for a new round

       current_round = pState.getRound();

       setUpRound(pState, pDue);

       std::cerr << "Round" << current_round << '\n';
       updateObservations(pState, observation_sequences);
     } else {

       updateObservations(pState, observation_sequences);

       setUpTimeStep(pState);
     }
     time_step++;


     if (time_step > TRAINING_START) {
       int most_certain_bird = -1;
       int movement_to_shoot = -1;
       double move_prob = 0;
       double highest_prob = SHOOT_PROB_THRESHOLD;

       for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
         if (pState.getBird(bird_index).isAlive()) {
           for (size_t movement = 0; movement < COUNT_MOVE; movement++) {

             move_prob = getMoveProb(bird_index, movement);

             if (move_prob > highest_prob) {
               if (guessSpecies(allAlphas[bird_index])[SPECIES_BLACK_STORK] < SHOOT_STORK_PROB_THRESHOLD ) {
                 most_certain_bird = bird_index;
                 movement_to_shoot = movement;
                 highest_prob = move_prob;
               }
             }
           }
         }
       }



       count_shooting_attempts++;

       std::cerr << "Time step: " << time_step << '\n';
       std::cerr << "Nr hit birds: " << count_bird_hits << '\n';
       if (most_certain_bird != -1) {
         std::cerr << "Shooting attempt nr: " << count_shooting_attempts << '\n';
         std::cerr << "Bird to shoot: " << most_certain_bird << '\n';
       }

       return Action(most_certain_bird, movement_map[movement_to_shoot]);

     } else {
       // This line choose not to shoot
       return cDontShoot;
     }

    //return cDontShoot;
    //This line would predict that bird 0 will move right and shoot at it
    // return Action(0, MOVE_RIGHT);
}

void scale(std::vector<double> & v) {
    double sum = 0;
    for (auto& val: v) {
      sum += val;
    }
    if (sum == 0) {
      return;
    }
    for (int index = 0; index < v.size(); index++) {
      v[index] = v[index] / sum;
    }
    return;
}

int getMin(const std::vector<double> & v){
  double min_val = 1000000000000000;
  int min_elem_index = -1;
  for (int col = 0; col < v.size(); col++) {
      if (v[col] < min_val){
         min_val = v[col];
         min_elem_index = col;
     }
    }
    return min_elem_index;
}

int getMax(const std::vector<double> & v){
  double max_val = -1000000000000000;
  int max_elem_index = -1;
  for (int col = 0; col < v.size(); col++) {
      if (v[col] > max_val){
         max_val = v[col];
         max_elem_index = col;
     }
    }
    return max_elem_index;
}

std::vector<double> guessSpecies(std::vector<std::vector<double>> lastAlphaVector) {

  // std::cerr << "Alpha: ";
  // for (auto prob1:  lastAlphaVector[0]) {
  //   std::cerr << prob1 << ", ";
  // }
  // std::cerr << '\n';

  std::vector<double> prob_dist_species = std::vector<double>(COUNT_SPECIES, 0);
  for (size_t species_index = 0; species_index < COUNT_SPECIES; species_index++) {

    double prob = observationProb(lastAlphaVector[species_index]);
    // std::cerr << "Prob: " << prob << '\n';
    prob_dist_species[species_index] = prob * prob_species_prior[species_index];

  }

  //
  //  std::cerr << "Before scale: ";
  //  for (auto prob1:  prob_dist_species) {
  //    std::cerr << prob1 << ", ";
  //  }
  //  std::cerr << '\n';
  // scale(prob_dist_species);
  // std::cerr << "After scale: ";
  // for (auto prob1:  prob_dist_species) {
  //   std::cerr << prob1 << ", ";
  // }
  // std::cerr << '\n';

  return prob_dist_species;
}

// bool haveAllSpecies(){
//   for (int species_index = 0; species_index < AllModels.size(); species_index++) {
//     if (AllModels[species_index].size() == 0) {
//       return false;
//     }
//   }
//   return true;
// }

// std::vector<int> getNotFoundSpecies(){
//   std::vector<int> notFoundSpecies;
//   for (size_t species = 0; species < AllModels.size(); species++) {
//     if (AllModels[species].size() == 0) {
//       notFoundSpecies.push_back(species);
//     }
//   }
//   return notFoundSpecies;
//
// }
//
// int confidentSpeciesGuess(std::vector<double> & guess_distribution){
//   double top_a = 1;
//   double top_b = 1;
//   int top_a_index = -1;
//   int top_b_index = -1;
//
//   for (size_t probability_index = 0; probability_index < guess_distribution.size(); probability_index++) {
//     if ((guess_distribution[probability_index] < top_a) || (guess_distribution[probability_index] < top_b)) {
//       if (top_a > top_b) {
//         top_a = guess_distribution[probability_index];
//         top_a_index = probability_index;
//       }else{
//         top_b = guess_distribution[probability_index];
//         top_b_index = probability_index;
//       }
//     }
//   }
//
//   if (top_a > top_b) {
//     double temp_a = top_a;
//     double temp_a_index = top_a_index;
//     top_a = top_b;
//     top_a_index = top_b_index;
//     top_b = temp_a;
//     top_b_index = temp_a_index;
//   }
//
//   if (top_a / top_b > GUESS_TOP_SIMILARITY_THRESHOLD) {
//     return SPECIES_UNKNOWN;
//   }
//   else{
//     return top_a_index;
//   }
//
// }

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */
     // bool have_all_species = haveAllSpecies();

     // std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);
     std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_PIGEON);
     if (current_round == 0) {
       return lGuesses;
     }
     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       std::vector<double> probabilities = prob_species[bird_index];

       double highest_prob = 0;
       int best_species_guess = -1;
       for (size_t species = 0; species < COUNT_SPECIES; species++) {
         // std::cerr << "Probability: " << probabilities[species] << '\n';
         if (probabilities[species] > highest_prob) {
           highest_prob = probabilities[species];
           best_species_guess = species;
         }

       }
       lGuesses[bird_index] = species_map[best_species_guess];
      }


     // for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
     //   std::vector<double> guess_distribution = guessSpecies(observation_sequences[(int) bird_index]);
     //   scale(guess_distribution);
     //
     //   // int best_guess_index = getMin(guess_distribution);
     //
     //
     //   int best_guess_index = confidentSpeciesGuess(guess_distribution);
     //   ESpecies guess = species_map[best_guess_index];
     //   if ((pState.getRound() == 0 || (!have_all_species && ((double) rand() / (RAND_MAX)) < PERCENTAGE_TO_GUESS_UNKNOWN)) && guess == SPECIES_UNKNOWN ) {
     //
     //     //just guess on species not revealed/found
     //     std::vector<int> notFoundSpecies = getNotFoundSpecies();
     //     int random_species_guess = notFoundSpecies[rand() % notFoundSpecies.size()];
     //
     //     guess = species_map[random_species_guess];
     //   }
     //   lGuesses[(int) bird_index] = guess;
     // }
     //
     // std::cerr << "Guesses: ";
     // for (auto  species:  lGuesses) {
     //   std::cerr << species << ", ";
     // }
     // std::cerr << '\n';

    return lGuesses;
}

void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */
    count_bird_hits += 1;
    // missed_birds[pBird]--;
    std::cerr << "HIT BIRD!!!" << std::endl;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */



     std::cerr << "Reveal: ";
     for (auto species:  pSpecies) {
       std::cerr << species << ", ";
     }
     std::cerr << '\n';



     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       int species_for_bird = pSpecies[bird_index];
       std::vector<int> & species_obs = observation_sequences_species[species_for_bird];
       std::vector<int> & bird_obs = observation_sequences[bird_index];
       species_obs.insert(species_obs.end(), bird_obs.begin(), bird_obs.end());

       species_observed[species_for_bird]++;
     }

     birds_observed += pState.getNumBirds();

     for (size_t species = 0; species < COUNT_SPECIES; species++) {
       std::vector<int> & species_obs = observation_sequences_species[species];
       if (species_obs.size() > 0) {
         std::cerr << "Training in progress..." << '\n';
         AllModels[species].train(species_obs);
         std::cerr << "Training Done" << '\n';

       }
       std::cerr << "Time left: " << pDue.remainingMs() << '\n';

       prob_species_prior[species] = (float) species_observed[species] / (float) birds_observed;
     }

     std::cerr << "Reveal done" << '\n';
     return;


}


} /*namespace ducks*/
