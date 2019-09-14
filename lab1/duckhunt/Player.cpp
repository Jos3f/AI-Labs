#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include "Matris.cp"
#include <vector>
#include <map>
#include <cmath>
// #include <r>

namespace ducks
{

// Forward declarations
struct HMM;
int guessSpecies(std::vector<int> observations);
int guessSpeciesAndBlackProb(std::vector<int> observations, double & black_stork_prob);


// Global mappings
std::map<int, EMovement> movement_map;
std::map<int, ESpecies> species_map;


std::vector<std::vector<int>> observation_sequences; // Bird, time


int time_step;
int current_round;

//best so far: score  with parameters
// int MAX_ITER_TRAIN = 40;
// int TRAINING_START = 70;
// double PROB_THRESHOLD = 0.8;
// double GUESS_LOGPROB_THRESHOLD = -5000;
// double PERCENTAGE_TO_GUESS_UNKNOWN = 0.2;
//

int MAX_ITER_TRAIN = 40;
int TRAINING_START = 70;
double PROB_THRESHOLD = 0.8;
double GUESS_LOGPROB_THRESHOLD = -5000;
double PERCENTAGE_TO_GUESS_UNKNOWN = 0.2;
double BLACK_STORK_LOGPROB_THRESHOLD = -4500;

std::vector<std::vector<HMM>> AllModels;



struct HMM {
  Matris a;
  Matris b;
  Matris pi;

  // int num_states;
  // int num_observations;

  int total_observations;

  std::vector<double> c;

  HMM(){
    a = Matris(5,5, {
      0.10, 0.20, 0.3, 0.15, 0.25,
      0.15, 0.10, 0.15, 0.35, 0.25,
      0.2, 0.08, 0.12, 0.25, 0.35,
      0.25, 0.35, 0.13, 0.07, 0.2,
      0.35, 0.15, 0.25, 0.12, 0.13
    });

    // a = Matris(3,3, {
    //   0.5, 0.20, 0.3,
    //   0.25, 0.40, 0.35,
    //   0.22, 0.52, 0.26
    // });

    b = Matris(9,5, {
      0.10, 0.12, 0.21, 0.15, 0.12, 0.05, 0.08, 0.09, 0.08,
      0.15, 0.10, 0.15, 0.08, 0.15, 0.1, 0.08, 0.12, 0.07,
      0.2, 0.08, 0.12, 0.05, 0.10, 0.05, 0.15, 0.2, 0.05,
      0.15, 0.15, 0.11, 0.07, 0.2, 0.08, 0.1, 0.02, 0.12,
      0.15, 0.15, 0.14, 0.12, 0.13, 0.1, 0.05, 0.05, 0.11
    });

    // b = Matris(9,3, {
    //   0.10, 0.12, 0.21, 0.15, 0.12, 0.05, 0.08, 0.09, 0.08,
    //   0.15, 0.10, 0.15, 0.08, 0.15, 0.1, 0.08, 0.12, 0.07,
    //   0.2, 0.08, 0.12, 0.05, 0.10, 0.05, 0.15, 0.2, 0.05
    // });

    pi = Matris(5,1, {
      0.12, 0.18, 0.2, 0.3, 0.2
    });

    // pi = Matris(3,1, {
    //   0.4, 0.32, 0.28
    // });

    // num_states = A.rows();
    // num_observations = B.cols();

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
  current_round = -1;
  // observation_sequences = vector<vector<int>>();
  // setUpRound();
  AllModels = std::vector<std::vector<HMM>>();
  for (int i = 0; i < COUNT_SPECIES; i++) {
    AllModels.push_back(std::vector<HMM>());
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
  observation_sequences = std::vector<std::vector<int>>();
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    observation_sequences.push_back(std::vector<int>());
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

std::vector<int> getDeathPenalty(const GameState & pState){
  int predMove = -1;
  int guiltyBird = -1;
  double highestProb = -1;

  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    if (pState.getBird(bird_index).isAlive()) {
      std::vector<int> observations = observation_sequences[(int) bird_index];
      double black_stork_prob;
      int likely_species = guessSpeciesAndBlackProb(observations, black_stork_prob);

      //conditions for death penalty
      if ((likely_species != SPECIES_BLACK_STORK)
      && (black_stork_prob < BLACK_STORK_LOGPROB_THRESHOLD)
      && (likely_species != SPECIES_UNKNOWN) ) {
        HMM birdModel = HMM();
        birdModel.train(observations);
        Matris lastAlpha = birdModel.getLastAlpha(observations);
        Matris movementDist = birdModel.calcNextObsDist(lastAlpha);

        for (int movement = 0; movement < movementDist.cols(); movement++) {
          if (movementDist(movement, 0) > highestProb && movementDist(movement, 0) > PROB_THRESHOLD) {
            highestProb = movementDist(movement, 0);
            guiltyBird = bird_index;
            predMove = movement;
          }
        }
      }
    }
  }
  return {guiltyBird, predMove};
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
     }
     time_step++;

     updateObservations(pState, observation_sequences);

     if (time_step > TRAINING_START) {
       std::vector<int> deathPenalty = getDeathPenalty(pState);
       int most_certain_bird = deathPenalty[0];
       int movement_to_shoot = deathPenalty[1];

       //if no birds next move is certain
       if (most_certain_bird == -1) {
         return cDontShoot;
       }

       std::cerr << "Time step: " << time_step << '\n';
       //std::cerr << "bird to shoot " << most_certain_bird << '\n';
       //std::cerr << "birds pred move " << movement_map[movement_to_shoot] << '\n';
       //std::cerr << "Current score " << pState.myScore() <<'\n';
       std::cerr << "Shooting attempt" << '\n';

       std::cerr << "Species to shoot " << guessSpecies(observation_sequences[(int) most_certain_bird]) << '\n';
       // int species_to_shoot
       // if (/* condition */) {
       //   /* code */
       // }

       return Action(most_certain_bird, movement_map[movement_to_shoot]);

     } else {
       // This line choose not to shoot
       return cDontShoot;
     }

    //return cDontShoot;
    //This line would predict that bird 0 will move right and shoot at it
    // return Action(0, MOVE_RIGHT);
}

int guessSpeciesAndBlackProb(std::vector<int> observations, double & black_stork_prob) {
  int best_guess = -1;
  double highest_similarity_prob = GUESS_LOGPROB_THRESHOLD;

  for (size_t speciesModels_index = 0; speciesModels_index < AllModels.size(); speciesModels_index++) {
    double avg_prob = -std::numeric_limits<double>::infinity();
    // std::cerr << "Start err prob" << '\n';
    for (size_t model_index = 0; model_index < AllModels[speciesModels_index].size(); model_index++) {
      if ((int) model_index == 0) {
        avg_prob = 0;
      }
      HMM hmm_model = AllModels[speciesModels_index][model_index];
      // std::cerr << "calcProb(observations): " << hmm_model.calcProb(observations) << '\n';
      double prob = hmm_model.calcProb(observations);
      // std::cerr << "Part avg_prob" << prob << '\n';

      avg_prob += prob / AllModels[speciesModels_index].size();
    }
    // std::cerr << "Avg_prob: " << avg_prob<< '\n';
    // std::cerr << "End err prob" << '\n';

    // std::cerr << "avg_prob: " <<  avg_prob << '\n';
    // std::cerr << "/* error message */" << '\n';
    if (avg_prob > highest_similarity_prob) {
      highest_similarity_prob = avg_prob;
      best_guess = speciesModels_index;
    }
    if ((int) speciesModels_index == (int) SPECIES_BLACK_STORK) {
      black_stork_prob = avg_prob;
    }
  }
  // std::cerr << "highest_similarity_prob: " << highest_similarity_prob << '\n';

  return best_guess;
}

int guessSpecies(std::vector<int> observations){

  double black_stork_trash;

  return guessSpeciesAndBlackProb(observations, black_stork_trash);
}

bool haveAllSpecies(){
  for (int species_index = 0; species_index < AllModels.size(); species_index++) {
    if (AllModels[species_index].size() == 0) {
      return false;
    }
  }
  return true;
}

std::vector<int> getNotFoundSpecies(){
  std::vector<int> notFoundSpecies;
  for (size_t species = 0; species < AllModels.size(); species++) {
    if (AllModels[species].size() == 0) {
      notFoundSpecies.push_back(species);
    }
  }
  return notFoundSpecies;

}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */
     bool have_all_species = haveAllSpecies();

     std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);
     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       ESpecies guess = species_map[guessSpecies(observation_sequences[(int) bird_index])];
       if ((pState.getRound() == 0 || (!have_all_species && ((double) rand() / (RAND_MAX)) < PERCENTAGE_TO_GUESS_UNKNOWN)) && guess == SPECIES_UNKNOWN ) {

         //just guess on species not revealed/found
         std::vector<int> notFoundSpecies = getNotFoundSpecies();
         int random_species_guess = notFoundSpecies[rand() % notFoundSpecies.size()];

         guess = species_map[random_species_guess];
       }
       lGuesses[(int) bird_index] = guess;
     }

     std::cerr << "Guesses: ";
     for (auto  species:  lGuesses) {
       std::cerr << species << ", ";
     }
     std::cerr << '\n';


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
     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       if (pSpecies[(int) bird_index] == -1) {
         continue;
       }
       HMM hmm_model;
       hmm_model.train(observation_sequences[(int) bird_index]);
       AllModels[pSpecies[(int) bird_index]].push_back(hmm_model);
     }
     std::cerr << "Reveal: ";
     for (auto  species:  pSpecies) {
       std::cerr << species << ", ";
     }
     std::cerr << '\n';

}


} /*namespace ducks*/
