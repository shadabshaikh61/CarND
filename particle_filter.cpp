/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	std::default_random_engine gen;
	num_particles = 10;
	is_initialized = true;
	
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);
	
	
	for(int i=0;i<num_particles;i++){
		Particle temp_particle;
		temp_particle.id = i;
		temp_particle.x = dist_x(gen);
		temp_particle.y = dist_y(gen);
		temp_particle.theta = dist_theta(gen);
		temp_particle.weight = 1;
		particles.push_back(temp_particle);
		
		weights.push_back(1);
	}
	return;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	if(yaw_rate!=0){
		for(int i=0;i<num_particles;i++){
			double x = particles[i].x;
			double y = particles[i].y;
			double theta = particles[i].theta;
			
			x += velocity * (sin(theta + yaw_rate * delta_t) - sin(theta))/yaw_rate;
			y += velocity * (cos(theta) - cos(theta + yaw_rate * delta_t))/yaw_rate;
			theta += yaw_rate*delta_t;
			
			particles[i].x = x;
			particles[i].y = y;
			particles[i].theta = theta;
		}
	}
	else{
		for(int i=0;i<num_particles;i++){
			double x = particles[i].x;
			double y = particles[i].y;
			double theta = particles[i].theta;
			
			x += velocity * (sin(theta));
			y += velocity * (cos(theta));
			
			particles[i].x = x;
			particles[i].y = y;
		}
	}
	return;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for(unsigned int i=0; i<observations.size(); i++){
		double x = observations[i].x;
		double y = observations[i].y;
		double temp = 999999; 
		LandmarkObs obs;
		for(unsigned int j=0; j<predicted.size(); j++){
			if( pow(predicted[j].x-x,2) + pow(predicted[j].y-y,2) < temp)
				obs = predicted[j];
				temp = pow(predicted[j].x-x,2) + pow(predicted[j].y-y,2)	;
		}
		observations[i].id = obs.id;		
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for(int i=0; i<num_particles; i++){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		double w = particles[i].weight;
		
		std::vector<LandmarkObs> transformed_lm;
		std::vector<LandmarkObs> predicted_inrange_lm;
		for(unsigned int j=0; j<observations.size(); j++){
			double xc = observations[j].x;
			double yc = observations[j].y;			

			LandmarkObs temp;
			temp.id = i*observations.size()+j;
			temp.x = cos(theta)*xc - sin(theta)*yc + x;
			temp.x = sin(theta)*xc + cos(theta)*yc + y;
			transformed_lm.push_back(temp);
		}
		
		for(unsigned int j=0; j<map_landmarks.landmark_list.size(); j++){
			double x_obs = map_landmarks.landmark_list[j].x_f;
			double y_obs = map_landmarks.landmark_list[j].y_f;
			if(pow(x_obs-x,2) + pow(y_obs-y,2) <= sensor_range){
				LandmarkObs temp_landmark;
				temp_landmark.x = x_obs;
				temp_landmark.y = y_obs;
				temp_landmark.id = map_landmarks.landmark_list[j].id_i;

				predicted_inrange_lm.push_back(temp_landmark);
			}
		}

		dataAssociation(predicted_inrange_lm,transformed_lm);
		
		double final_weight = 1;
		for(unsigned int k=0; k<observations.size(); k++){
			LandmarkObs temp;
			for(unsigned int j=0;j<predicted_inrange_lm.size();j++){
				if(predicted_inrange_lm[j].id == observations[k].id){
					temp = observations[k];
					break;
				}
			}
			double x_obs = observations[k].x;
			double y_obs = observations[k].y;
			final_weight*=exp(-1*(pow(x_obs-temp.x,2)/2*pow(std_landmark[0],2)+pow(y_obs-temp.y,2)/2*pow(std_landmark[1],2)));
			final_weight/=2*M_PI*std_landmark[0]*std_landmark[1];

		}
		//std::cout<<"Particle : "<<i<<"\tweight :"<<w<<endl;
		particles[i].weight = final_weight;
		weights[i] = final_weight;
	}	
/*
	double sum_w=0;
	for(int i=0;i<num_particles;i++){
		sum_w+=particles[i].weight;
	}

	for(int i=0;i<num_particles;i++){
		particles[i].weight/=sum_w;
		weights[i]/=sum_w;
	}*/

	return;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	std::vector<Particle> new_particle_list;	

	std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(),weights.end());
  for(int n=0; n<num_particles; ++n) {
		int i = d(gen);
		new_particle_list.push_back(particles[i]);
  }
	std::vector<Particle> particles(new_particle_list);
	return;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
		return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
