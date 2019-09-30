#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include "Matris.cp"
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <time.h>

/**
 * Algorithm explanation
 *
 * We start by not shooting anything the first round
 *
 * After the first round, we guess that each bird is species pigeon
 *
 * After each reveal, we take the observations for the birds in that round. These new observations are appended to
 * the stored observations for that species and then the species model is re-trained. Each species has one stored model for an environment.
 *
 * For shooting: For each bird and the 9 movements, we calculate how probable it is to do that move. This is done by utilizing the 6 species models in the following way:
 * 1. calculate the probability of seeing that move, given all previous observations and the species model.
 * 2. multiply the result in (1.) with the probability of the species for the bird
 * 3. Pick the bird that makes the move we are most certain about, given that the probability for it being the stork is
 * low enough and that the probability of that move is larger than our threshold for shooting.
 *
 * Guesses in future rounds are done by calculating the probability of observation sequences for each species model and then multiply by the prior probability of that species.
 * The species corresponding to the highest such probability becomes the guess. Note that we always guess on some species here.
 *
 *
 */

namespace ducks
{

// Forward declarations
struct HMM;
std::vector<double> guessSpecies(std::vector<std::vector<double>> lastAlphas);
void scale(std::vector<double> & v);
void updateObservations(const GameState &pState, std::vector<std::vector<int>> &observation_sequences);

// Global mappings, used for converting int to EMovement/ESpecies only
std::map<int, EMovement> movement_map;
std::map<int, ESpecies> species_map;

std::vector<std::vector<int>> observation_sequences; // Bird, time. All observations so far in this round
std::vector<std::vector<int>> observation_sequences_species; // species, time. All observations so far in this Environment for each species
std::vector<int> species_observed; // Counter for how many times a species has occured so far in this environment
int birds_observed;                 // Total birds observed so far in this environment
std::vector<double> prob_species_prior; // prior probabilities for each species

std::vector<std::vector<std::vector<double>>> allAlphas; // birds in round, species, pattern/hidden. Alphas for each bird this round and for each model in the species model matrix (AllModels).
std::vector<std::vector<double>> prob_species;          // Birds in round, species. Probability distribution for what species each bird is, each round.

int time_step;
int current_round;
int count_bird_hits = 0; // For debugging
int count_shooting_attempts = 0; // For debugging


int MAX_ITER_TRAIN = 100;
int TRAINING_START = 70;
double SHOOT_PROB_THRESHOLD = 0.5; // How certain we have to be that a move is in a specific direction
double SHOOT_STORK_PROB_THRESHOLD = 0.3;  // How certain we have to be that the bird is a stork to not shoot
int NUM_STATES = 4; // Hidden states in all HMMs

std::vector<HMM> AllModels; // All models so far, max one per species

struct HMM {
  Matris a;
  Matris b;
  Matris pi;

  std::vector<double> c; // scaling

  HMM(int num_states){

    a = generate_a_matris(num_states);

    b = generate_b_matris(num_states);

    pi = generate_pi(num_states);
  }

  /**
   * Random value for element in a matrix
   */
  double generate_element_val(int num_cols){

      double uni_rand = (double)rand() / RAND_MAX;
      double scale_factor =  (double) num_cols;
      double avg_element_value = 1/scale_factor;
      double min_val = (avg_element_value - avg_element_value/2);
      double max_val = (avg_element_value + avg_element_value/2);

      double element_val = min_val + uni_rand * (max_val - min_val);
      return element_val;
  }

  /**
   * Matrix with random values, each row sums to 1
   * @param num_cols
   * @param num_states
   * @return
   */
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

  /**
   * Observation probability distribution for the next step given current state probability distribution
   * @param current_pi
   * @return
   *
  Matris calcNextObsDist(Matris current_pi){
    return current_pi * a * b;
  }
*/

  /**
   * Calculate the aplha matrix for the given observation
   * @param observations
   * @return
   */
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

    // Avoid NaN and inf in matrices, quick hack
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

      // Avoid NaN and inf in matrices, quick hack
      if (isinf(c[t]) || isnan(c[t]) ) {
        c[t] = 10e+50;
      }

      for (int i = 0; i < a.rows(); i++) {
        allAlpha(t, i) *= c[t];
      }
    }

    return allAlpha;
  }

  /**
   * Calculate all beta values for the observation
   * @param observations
   * @return
   */
  Matris calcBeta(const std::vector<int> &observations){
    int total_observations = observations.size();
    Matris allBeta(total_observations, a.rows());
    // calculate Beta pass T
    for (int state = 0; state < a.rows(); state++) {
      allBeta(total_observations - 1, state) = 1 * c[total_observations - 1];
    }

    for (int time_step = total_observations - 2; time_step >= 0; time_step--) {
      for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.rows(); j++) {
          allBeta(time_step, i) += a(j, i) * b(observations[time_step + 1], j) * allBeta(time_step + 1, j);
        }
        // Scale beta t with same factor as alpha t
        allBeta(time_step, i) *= c[time_step];
      }
    }
    return allBeta;
  }

  /**
   * Smoothens row in matrix by making almost zero values a bit larger, prevents NaN and inf in future calculations
   * @param row
   * @param matris
   */
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

   /**
    * Train the HMM using baum-welch algorithm with the provided observations
    * @param observations
    */
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
            a(j,i) = numer/denom;
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
            b(k,j) = numer/denom;
        }
      }

      // Compute log[P(O|λ)] to use as stop training criteria
      oldLogProb = logProb;
      logProb = 0;
      for (int time_step = 0; time_step < total_observations; time_step++) {
        logProb -= log(c[time_step]);
      }
      itertation++;
    }

    for (int row = 0; row < a.rows(); row++) {
      smoothRow(row, a);
    }

    for (int row = 0; row < b.rows(); row++) {
      smoothRow(row, b);
    }

    for (int row = 0; row < pi.rows(); row++) {
      smoothRow(row, pi);
    }

    return;
  }

  /**
   * Returns the sum of the alphas at the last time step
   * @param allAlpha
   * @return
   */
  double sumAlphaT(const Matris allAlpha){
    double sum = 0;
    for (int state = 0; state < a.rows(); state++)
      sum += allAlpha(allAlpha.cols() - 1, state);
    return sum;
  }

};

/**
 * This is called each time we start a new environment
 */
Player::Player()
{
  // srand(time(NULL));

  // Initialize variables needed for the environment
  current_round = -1;
  species_observed = std::vector<int>(COUNT_SPECIES, 0);
  prob_species_prior = std::vector<double>(COUNT_SPECIES, 0);

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

/**
 * Do everything that is needed for a new round
 * @param pState
 * @param pDue
 */
void setUpRound(const GameState &pState, const Deadline &pDue){

    // Initialize variables needed for each round
  time_step = -1;
  prob_species = std::vector<std::vector<double>>(pState.getNumBirds(), std::vector<double>(COUNT_SPECIES, 0));

  // Calculate the alphas and probabilities for species for each bird the very first time step
  allAlphas = std::vector<std::vector<std::vector<double>>>(pState.getNumBirds(), std::vector<std::vector<double>>(COUNT_SPECIES, std::vector<double>(NUM_STATES, 0)));
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    for (size_t species = 0; species < COUNT_SPECIES; species++) {
      for (size_t state = 0; state < NUM_STATES; state++) {
        allAlphas[bird_index][species][state] = AllModels[species].b(pState.getBird(bird_index).getLastObservation(), state) * AllModels[species].pi(state, 0);
      }
    }
    prob_species[bird_index] = guessSpecies(allAlphas[bird_index]);
  }

  observation_sequences = std::vector<std::vector<int>>();
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    observation_sequences.push_back(std::vector<int>());
  }
  updateObservations(pState, observation_sequences);
}

/**
 * Updates the observations matrix for each bird using the latest observations, if the bird is alive. Should be called each time step
 * @param pState
 * @param observation_sequences
*/
void updateObservations(const GameState &pState, std::vector<std::vector<int>> &observation_sequences){
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    Bird bird = pState.getBird((int) bird_index);

    if (bird.isAlive()) {
      observation_sequences[(int) bird_index].push_back(bird.getLastObservation());
    }
  }
}

/**
 * Calculates the alphas in the next step using a model, the previous alphas and the last observation
 * @param model -
 * @param alphas -
 * @param observation -
 * @return
 */
std::vector<double> alphasNextStep(HMM model, std::vector<double> alphas, int observation){
  std::vector<double> alphas_next_step(NUM_STATES, 0);
  for (size_t state_from = 0; state_from < NUM_STATES; state_from++) {
    for (size_t state_to = 0; state_to < NUM_STATES; state_to++) {
      alphas_next_step[state_to] += alphas[state_from] * model.a(state_to, state_from)*model.b(observation, state_to);
    }
  }
  return alphas_next_step;
}

/**
 * Updates our global allAlphas matrix using alphasNextStep for each bird and species
 * @param pState
 */
void nextAlphaUpdate(GameState pState){
  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    if (pState.getBird(bird_index).isAlive()) {
      for (size_t species = 0; species < COUNT_SPECIES; species++) {
        allAlphas[bird_index][species] = alphasNextStep(AllModels[species], allAlphas[bird_index][species], observation_sequences[bird_index][observation_sequences[bird_index].size() - 1]);
      }
    }
  }
}

/**
 * This is called each time step
 * @param pState
 */
void setUpTimeStep(GameState pState){
  nextAlphaUpdate(pState);

  for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
    prob_species[bird_index] = guessSpecies(allAlphas[bird_index]);
  }
  return;
}

/**
 * Computes probability of observations, given the model
 * @param lastAlphaVector
 * @return
 */
double observationProb(std::vector<double> lastAlphaVector){
  double prob = 0;
  for (size_t state = 0; state < NUM_STATES; state++) {
    prob += lastAlphaVector[state];
  }
  return prob;
}

/**
 * Computes probability of a certain movement for a bird
 * @param bird_index
 * @param movement
 * @return
 */
double getMoveProb(int bird_index, int movement){
  double probability = 0.000;
  for (size_t species = 0; species < COUNT_SPECIES; species++) {
    // First computes alphas one step ahead and then uses it in observationProb
    double prob_next = observationProb(alphasNextStep(AllModels[species], allAlphas[bird_index][species], movement));
      // First computes alphas this time step and then uses it in observationProb to get prob of all observations given model
    double prob_this_time_step = observationProb(allAlphas[bird_index][species]);
    if (prob_this_time_step > 0) {
      //The quotient gives probability of next observation (movement) given model (includes species guess) and then multiplies by prob of species to take in that information
      probability += ((prob_next / prob_this_time_step) * prob_species[bird_index][species]);
    }
  }
  return probability;
}

/**
  * Here you should write your clever algorithms to get the best action.
  * This skeleton never shoots.
  */
Action Player::shoot(const GameState &pState, const Deadline &pDue)
{

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

     // Don't want to shoot too early when there is not much data
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
               // better safe than sorry
               if (prob_species[bird_index][SPECIES_BLACK_STORK] < SHOOT_STORK_PROB_THRESHOLD ) {
                 most_certain_bird = bird_index;
                 movement_to_shoot = movement;
                 highest_prob = move_prob;
               }
             }
           }
         }
       }




       std::cerr << "Time step: " << time_step << '\n';
       std::cerr << "Nr hit birds: " << count_bird_hits << '\n';
       if (most_certain_bird != -1) {
           count_shooting_attempts++;
           std::cerr << "Shooting attempt nr: " << count_shooting_attempts << '\n';
         std::cerr << "Bird to shoot: " << most_certain_bird << '\n';
       }

       return Action(most_certain_bird, movement_map[movement_to_shoot]);

     } else { // If not above training start
       // This line choose not to shoot
       return cDontShoot;
     }

}

/**
 * Scale array s.t. the row sums to one
 * @param v
 */
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

/**
 * Returns probability distribution of all species for a bird given its alpha matrix. Includes information about prior probability of species
 * @param lastAlphaVector
 * @return
 */

std::vector<double> guessSpecies(std::vector<std::vector<double>> lastAlphaVector) {

  std::vector<double> prob_dist_species = std::vector<double>(COUNT_SPECIES, 0);
  for (size_t species_index = 0; species_index < COUNT_SPECIES; species_index++) {

    double prob = observationProb(lastAlphaVector[species_index]);
    prob_dist_species[species_index] = prob * prob_species_prior[species_index];

  }
  // scale to make probs comparable
  scale(prob_dist_species);
  return prob_dist_species;
}

/**
 * Here you should write your clever algorithms to guess the species of each bird.
 * This skeleton makes no guesses, better safe than sorry!
 */
std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{

    // Default guess PiGEON worked well
     std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_PIGEON);
     if (current_round == 0) { // cant make educated guesses here
       return lGuesses;
     }

     //Guess on species with highest probability, ALWAYS guess
     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       std::vector<double> probabilities = prob_species[bird_index];

       double highest_prob = 0;
       int best_species_guess = -1;
       for (size_t species = 0; species < COUNT_SPECIES; species++) {
         if (probabilities[species] > highest_prob) {
           highest_prob = probabilities[species];
           best_species_guess = species;
         }

       }
       lGuesses[bird_index] = species_map[best_species_guess];
      }

    return lGuesses;
}
/**
 * If you hit the bird you are trying to shoot, you will be notified through this function.
 */
void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    count_bird_hits += 1;
    std::cerr << "HIT BIRD!!!" << std::endl;
}

/**
 * If you made any guesses, you will find out the true species of those birds in this function.
 *
 */
void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{

    //Debugging
     std::cerr << "Reveal: ";
     for (auto species:  pSpecies) {
       std::cerr << species << ", ";
     }
     std::cerr << '\n';


     for (size_t bird_index = 0; bird_index < pState.getNumBirds(); bird_index++) {
       int species_for_bird = pSpecies[bird_index];
       std::vector<int> & species_obs = observation_sequences_species[species_for_bird]; //creates reference to stored observations for the species of the bird (environment unique)
       std::vector<int> & bird_obs = observation_sequences[bird_index]; //creates reference to observations for the species of the bird (round unique)
       species_obs.insert(species_obs.end(), bird_obs.begin(), bird_obs.end()); // concat stored observation sequences with new ones

       species_observed[species_for_bird]++;
     }

     birds_observed += pState.getNumBirds();

     // Re-train the stored models with all species observations from all rounds up until now
     for (size_t species = 0; species < COUNT_SPECIES; species++) {
       std::vector<int> & species_obs = observation_sequences_species[species];
       //Train only if the species has been observed at least once, even doe its only needed when we observe a species in this new round
        if (species_obs.size() > 0) {
         std::cerr << "Training in progress..." << '\n';
         AllModels[species].train(species_obs);
         std::cerr << "Training Done" << '\n';

       }
       std::cerr << "Time left: " << pDue.remainingMs() << '\n';

        // update species prior as the empirical distribution
       prob_species_prior[species] = (float) species_observed[species] / (float) birds_observed;
     }

     std::cerr << "Reveal done" << '\n';
     return;


}


} /*namespace ducks*/
